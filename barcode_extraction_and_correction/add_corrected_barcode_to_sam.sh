#!/bin/bash
#
# Purpose:
#   Append corrected barcode SAM tags ("CR:Z:", "CY:Z:", "CB:Z:") from a corrected
#   barcode file to each record in a SAM file.
#
#   The corrected barcode file is expected to have read names matching those in
#   the SAM file. Multiple SAM records with the same read name may be adjacent;
#   once a new read name is seen, the previous one will never appear again.


set -e
set -o pipefail


decompress_cat_cmd='cat';
decompress_gzip_cmd='gzip -c -d';
decompress_igzip_cmd='igzip -c -d';
decompress_zstd_cmd='zstd -c -d';


add_corrected_barcode_to_sam () {
    local sam_filename="${1}";
    local corrected_bc_filename="${2}";
    local sam_output_filename="${3}";

    if [ ${#@} -lt 3 ] ; then
        printf '\nUsage:\n';
        printf '    add_corrected_barcode_to_sam \\\n';
        printf '        sam_input \\\n';
        printf '        corrected_bc_filename \\\n';
        printf '        sam_output \\\n\n';
        printf 'Purpose: Append raw barcode, raw barcode quality and corrected barcode SAM tags\n';
        printf '         ("CR:Z:", "CY:Z:", "CB:Z:") to each record in a SAM file.\n\n';
        printf 'Parameters:\n';
        printf '  - sam_input:             Input SAM filename ("-" or "stdin" for standard input).\n';
        printf '  - corrected_bc_filename: File with corrected barcodes (one per line, FASTQ read\n';
        printf '                           name format with SAM tags as comment).\n';
        printf '                           Can be gzip (.gz) or zstd (.zst) compressed.\n';
        printf '  - sam_output:            Output SAM filename ("-" or "stdout" for standard output).\n\n';
        return 1;
    fi


    if type igzip > /dev/null 2>&1 ; then
        local decompress_gzipped_cmd="${decompress_igzip_cmd}";
    else
        local decompress_gzipped_cmd="${decompress_gzip_cmd}";
    fi

    # Detect if input SAM file is compressed or stdin.
    if [ "${sam_filename}" = "-" ] || [ "${sam_filename}" = "stdin" ] ; then
        local read_sam_cmd="cat /dev/stdin";
    elif [ "${sam_filename}" != "${sam_filename%.gz}" ] ; then
        local read_sam_cmd="${decompress_gzipped_cmd} ${sam_filename}";
    elif [ "${sam_filename}" != "${sam_filename%.zst}" ] ; then
        local read_sam_cmd="${decompress_zstd_cmd} ${sam_filename}";
    else
        local read_sam_cmd="${decompress_cat_cmd} ${sam_filename}";
    fi

    # Detect if input corrected barcode file is compressed (gzip/zstd) or not.
    if [ "${corrected_bc_filename}" != "${corrected_bc_filename%.gz}" ] ; then
        local decompress_corrected_bc_file_cmd="${decompress_gzipped_cmd}";
    elif [ "${corrected_bc_filename}" != "${corrected_bc_filename%.zst}" ] ; then
        local decompress_corrected_bc_file_cmd="${decompress_zstd_cmd}";
    else
        local decompress_corrected_bc_file_cmd="${decompress_cat_cmd}";
    fi

    # Detect output destination.
    if [ "${sam_output_filename}" = "-" ] || [ "${sam_output_filename}" = "stdout" ] ; then
        local write_sam_cmd="cat > /dev/stdout";
    else
        local write_sam_cmd="cat > ${sam_output_filename}";
    fi


    if ! type mawk > /dev/null 2>&1 ; then
        printf 'Error: "mawk" not found or executable.\n';
        return 1;
    fi

    if ! type "${decompress_corrected_bc_file_cmd%% *}" > /dev/null 2>&1 ; then
        printf 'Error: "%s" not found or executable.\n' "${decompress_corrected_bc_file_cmd%% *}";
        return 1;
    fi


    mawk \
        -v sam_filename="${sam_filename}" \
        -v corrected_bc_filename="${corrected_bc_filename}" \
        -v sam_output_filename="${sam_output_filename}" \
        -v read_sam_cmd="${read_sam_cmd}" \
        -v decompress_corrected_bc_file_cmd="${decompress_corrected_bc_file_cmd}" \
        -v write_sam_cmd="${write_sam_cmd}" \
    '
    BEGIN {
        read_corrected_bc_file_cmd = decompress_corrected_bc_file_cmd " " corrected_bc_filename;

        read_name_corrected_bc = "";
        corrected_bc_sam_tags = "";
        corrected_bc_line = "";

        # Read SAM file line by line.
        while ( (read_sam_cmd | getline sam_line) > 0 ) {

            # Pass through SAM header lines unchanged.
            if ( substr(sam_line, 1, 1) == "@" ) {
                print sam_line | write_sam_cmd;
                continue;
            }

            # Extract read name from SAM record (first tab-delimited field).
            tab_pos = index(sam_line, "\t");
            if (tab_pos > 1) {
                read_name_sam = substr(sam_line, 1, tab_pos - 1);
                sam_rest = substr(sam_line, tab_pos);
            } else {
                print "Error: Invalid SAM record: \"" sam_line "\"." > "/dev/stderr";
                exit(1);
            }

            # Remove CB and UMI suffix from read name, if present.
            sub(/:[ACGTN+]+$/, "", read_name_sam);

            # Advance the corrected barcode file until we find a matching read name.
            # Because both files are in the same read order, and multiple SAM records
            # can share the same read name, we only advance the barcode file when the
            # current barcode read name does not match the SAM read name.
            while ( read_name_corrected_bc != read_name_sam ) {
                if ( (read_corrected_bc_file_cmd | getline corrected_bc_line) > 0 ) {
                    # Split corrected barcode line on spaces:
                    #   - read name: strip "@" from the start.
                    #   - corrected barcode info as SAM tags separated by "\t":
                    #       CR:Z:raw_bc\tCY:Z:raw_bc_qual\tCB:Z:corrected_bc_seq

                    # Find first space position in the corrected barcode line.
                    read_name_corrected_bc_space_pos = index(corrected_bc_line, " ");

                    if (read_name_corrected_bc_space_pos > 0) {
                        read_name_corrected_bc = substr(corrected_bc_line, 2, read_name_corrected_bc_space_pos - 2);
                        bc_sam_tags = substr(corrected_bc_line, read_name_corrected_bc_space_pos + 1);
                    } else {
                        read_name_corrected_bc = substr(corrected_bc_line, 2);
                        bc_sam_tags = "";
                    }

                    # Remove CB and UMI suffix from read name, if present.
                    sub(/:[ACGTN+]+$/, "", read_name_corrected_bc);

                } else {
                    # Corrected barcode file is exhausted but SAM records remain.
                    print "Error: Corrected barcode file ran out of records before SAM file (SAM read name: \"" read_name_sam "\")." > "/dev/stderr";
                    exit(1);
                }
            }

            # Append corrected barcode SAM tags (tab-separated) to the SAM record.
            if ( bc_sam_tags != "" ) {
                # Append SAM tags a the end of the SAM line (assume that SAM tags are unique in the whole line).
                print sam_line "\t" bc_sam_tags | write_sam_cmd;
            } else {
                print sam_line | write_sam_cmd;
            }
        }

        # Close open file handles.
        close(read_sam_cmd);
        close(read_corrected_bc_file_cmd);
        close(write_sam_cmd);
    }'

    return $?
}



add_corrected_barcode_to_sam "${@}";
