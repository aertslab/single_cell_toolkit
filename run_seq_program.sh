#!/bin/bash

# Set Seq root dir from environment variable (else use hard coded path).
seq_root_dir="${seq_root_dir:-/staging/leuven/stg_00002/lcb/ghuls/software/seq}"


run_seq_program () {
    if [ ${#@} -lt 1 ] ; then
        printf 'Usage: run_seq_program program.seq [program_args]\n';
        return 1;
    fi

    OMP_NUM_THREADS=4 LIBRARY_PATH="${seq_root_dir}/lib/seq:${LIBRARY_PATH}" "${seq_root_dir}/bin/seqc" run -release "${@}";
}



run_seq_program "${@}";
