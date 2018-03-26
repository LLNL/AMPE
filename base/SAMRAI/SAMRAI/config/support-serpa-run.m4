dnl Define a macro for supporting generalized serial-parallel run.


AC_DEFUN([CASC_SUPPORT_SERPA],[
dnl Support a generalized way to run a program in serial or parallel mode.
dnl (serpa is serial + parallel).
dnl
dnl Arg1: name (may include relative path) to give to shell script created.
dnl   If omitted, the name will be serpa-run, in the top level directory.
dnl
dnl This macro brings out an automake bug.  See "IMPORTANT NOTE:" below
dnl for comments and how to circumvent.

# Begin macro $0


# The generalized script needs to know about the system
# so it can form the appropriate parallel run command line.
AC_REQUIRE([AC_CANONICAL_TARGET])


dnl The parallel run program (usually mpirun) is hardwired for now,
dnl but it should be more flexible.
dnl I plan to eventualy use the following variable,
dnl but it is not currently used.
AC_ARG_WITH(parallel-run-bin,[
  --with-parallel-run-bin=STRING
		Specify the parallel run binary (i.e., mpirun, dmpirun,
		poe, etc) to be used in the serial-parallel run script.
],[
case "$with_parallel_run_bin" in
  no)
    # Set PARALLEL_RUN_BIN to blank to remove it from the serpa-run script.
    PARALLEL_RUN_BIN=;;
  yes)
    # Unset the variable PARALLEL_RUN_BIN so it is guessed in the next step.
    unset PARALLEL_RUN_BIN;
    ;;
  *)
    # Set PARALLEL_RUN_BIN to user specification.
    PARALLEL_RUN_BIN="$with_parallel_run_bin"
    ;;
esac
CASC_AC_LOG_VAR(with_parallel_run_bin PARALLEL_RUN_BIN target_os, with-parallel-run-bin given)

],[
unset PARALLEL_RUN_BIN;
CASC_AC_LOG_VAR(with_parallel_run_bin PARALLEL_RUN_BIN target_os, with-parallel-run-bin NOT given)
])
# If PARALLEL_RUN_BIN is unset, guess it.
if test ! "${PARALLEL_RUN_BIN+set}" = set; then
  case "$target_os" in
    osf*) PARALLEL_RUN_BIN=dmpirun
    ;;
    *) PARALLEL_RUN_BIN=mpirun
    ;;
  esac
  CASC_AC_LOG_VAR(with_parallel_run_bin PARALLEL_RUN_BIN target_os, after setting PARALLEL_RUN_BIN)
fi
AC_SUBST(PARALLEL_RUN_BIN)


define(btng_serpa_run_fn,ifelse($1,,serpa-run,$1))
dnl Create serpa-run.in file (at autoconf time, before configure time)
dnl File serpa-run.in will be used at configure time to create serpa-run
dnl in the compile tree.
syscmd([ cat <<'__EOM__'> ]btng_serpa_run_fn[.in
#!/bin/sh

# This script runs an executable program serially or in parallel.

# The syntax for running programs in parallel differs from platform to platform,
# software to software.  This script automatically uses the correct syntax for
# the environments that it knows about.  For the environments it does not know,
# it has a default that may work.

# This script defines the set of tests to put the program through.
# This set is parametrized by these variables:
# nproc_list: space- or comma-separated list of number of processors to use.


# The name and directory of this script
# (Do not rely on existence of dirname and basename programs.)
script_name=`echo ${0} | sed -e 's:.*/::'`;
dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;

# Define a way to gracefully die from within this script.
# The die function is similar to Perl's.
# Arguments are: <exit_value> <exit_message>
die () {
    echo "ERROR_MESSAGE_FROM ${script_name}:"
    if [ -n "${2}" ]; then echo "Error ${2}"; fi
    if [ -n "${1}" ]; then exit ${1}; fi
    exit 99
}


# When no argument is given, print the help message and exit.
if [ "${1}" = '-h' ] || [ "${1}" = '--help' ] || [ ${#} -lt 2 ] ; then
  cat <<-_EOF_
	Usage: ${script_name} <list of nproc> <program name>

	This is a generalized script to run a program in serial
	and/or parallel.

	The number of processors are given by <list of nproc>,
	which must be a comma- or space-delimited list of integers.
	Example: "${script_name} 1,2,5 parallel_program" runs
	parallel_program 3 times with 1, 2 and 5 processors in turn.

	The special case of number of processors = 0 means to run serially.
	This differs from one processor in that one processor means to run
	in parallel, but with one processor.

	This script exits with an error if any instance of
	running the program fails.

	If SERPA_REDIRECT_OUTPUT_TO is defined in the environment,
	standard output of the program is redirected there.
	If SERPA_REDIRECT_ERRORS_TO is defined in the environment,
	error output of the program is redirected there.
	Outputs directly from THIS script are NOT affected
	by these environment variables.

	If SERPA_MAX_FAILS is defined, ${script_name} will
	continue until that many failures before exiting.
	The default is 1 failure (${script_name} exits on the
	first failure).  The exit value will always be the
	failure count.

	When a serial run os made, if SERPA_PERFORMANCE_FILE
	is defined and is the name of a file, that file is searched
	for performance data.  The run is timed and the result
	is compared to the performance data found.
	_EOF_
exit
fi

# How to invoke a GNU compatible time command
if [ -z "${SERPA_GNUTIME}" ]; then
    SERPA_GNUTIME=/usr/bin/time
fi

# Percent difference for performance comparisons
if [ -z "${SERPA_DIFFERENCE}" ]; then
    SERPA_DIFFERENCE=0.1
fi


# Check the redirection environment variables.
# See the help message for how these are used.
if [ -n "${SERPA_REDIRECT_OUTPUT_TO}" ]; then
  output_redirection_string="1> ${SERPA_REDIRECT_OUTPUT_TO}"
else
  unset output_redirection_string
fi
if [ -n "${SERPA_REDIRECT_ERRORS_TO}" ]; then
  errors_redirection_string="2> ${SERPA_REDIRECT_ERRORS_TO}"
else
  unset errors_redirection_string
fi

# Save the arguments for later reference.
serpa_args="${@}"


# Get the number of processors from the first argument of this script.
nproc_list=${1}; shift;
# We allow comma-delimited lists, so we now remove those commas.
echo ${nproc_list} | grep '^[0-9 ,]\{1,\}$' > /dev/null	\
  || die 1 "Invalid number of processors list '${nproc_list}'"
nproc_list=`echo ${nproc_list} | sed 's/,/ /g'`

# Get the program name from the next argument of this script.
program=${1}; shift;
test -x ${program} || die 1 "No execute permission on '${program}'"


# Determine host name for use below.
if [ ${?}{HOST} ]; then
  HOST=`uname -n`
  export HOST
fi


# Use additional_env to specify additional environments that
# should be set before running.  Although you can, do not use
# this to set environments for programs.  It is meant to be
# environments for the parallel execution program such as
# mpirun.
additional_env=


# Variables required for using parallelrun.  These may not be needed,
# but if they are needed and unset, there will be an error.
# Some of these strings are determined at configure time.
parallelrun_prog="@PARALLEL_RUN_BIN@";	# Name of parallelrun program.


# Determine what platform we are on.
target_cpu=@target_cpu@
target_os=@target_os@
target_vendor=@target_vendor@


# Determine the machine file name.
serpa_machine_file='serpa.machines'


# Define a function to run a program in parallel.
# We have to do this because different environments require different syntaxes.
# The function arguments are:
#   program name
#   number of processor
#
# We define the function in a big if-else statement based on the
# target computer.  We assume that there is a one-to-one correspondence.
# If there is not, there may have to be a nested if-structure.
#
if echo "${HOST}" | egrep '^(blue|frost)[0-9]' > /dev/null ; then

  # For blue, a specific singleton platform.
  run_multiproc () {
    MP_RESD="YES"
    MP_HOSTFILE=""
    MP_EUILIB=us
    MP_EUIDEVICE=css0
    export MP_RESD MP_HOSTFILE MP_EUILIB MP_EUIDEVICE
  # IBM shell functions eat their parameters after the first assignment
  # from them so save the parameters first.
    case "${HOST}" in
      blue*) proc_per_node=4;;
      frost*) proc_per_node=16;;
    esac
    program=${1}; nproc=${2}; shift 2
    nodes=`expr 1 + \( ${nproc} - 1 \) / ${proc_per_node}`
    com="${additional_env} ${parallelrun_prog} ${program}"
    com="${com} -rmpool 0 -nodes ${nodes} -procs ${nproc} ${@}"
    # com="${program} ${@} -procs ${nproc}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${HOST}" | grep '^tc2k' > /dev/null ; then

  # For tc2k, a specific singleton platform (contributed by Brian Miller).
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} prun -n ${nproc} ${program} ${@}"
    com="${program} ${@} -p ${nproc}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${HOST}" | egrep '^(mcr|pengra|thunder|zeus)[0-9]' > /dev/null ; then

  # For LC Linux clusters using srun.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} srun -n${nproc} -N${nproc} -p pdebug ${program} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${target_os}" | grep '^osf'	> /dev/null ; then

  # For Dec OSF.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} dmpirun -np ${nproc} ${program} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

elif echo "${target_os}" | grep '^solaris' > /dev/null	\
  || echo "${target_os}" | grep '^irix' > /dev/null	\
  || echo "${target_os}" | grep '^linux' > /dev/null	\
  ; then

  # Most platforms fall into this case.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} ${parallelrun_prog}"
    com="${com} -machinefile ${dir_name}/${serpa_machine_file} -np ${nproc} ${program} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

else

  # Simplest case.
  # This is generic.  It may not work, but it is our best guess
  # without knowledge of the system.
  run_multiproc () {
    program=${1}; nproc=${2}; shift 2
    com="${additional_env} ${program} -np ${nproc} ${@}"
    echo ${com}
    eval ${com} ${output_redirection_string} ${errors_redirection_string}
    return ${?}
  }

fi


# Initialize the failure count.
serpa_num_failures=0
test -z "${SERPA_MAX_FAILS}" && SERPA_MAX_FAILS=1

# Run the program.
for nproc in ${nproc_list}; do
  cat <<-_EOM_
	${script_name}================================================
	RUNNING: ${script_name} ${nproc} ${program} ${@}
	${script_name}::::::::::::::::::::::::::::::::::::::::::::::::
	_EOM_
  if test "${nproc}" = 0 ; then
    # Run serially.

    # If performance file exists then with performance testing
    # otherwise do a normal sequential run
    if test "${SERPA_PERFORMANCE_FILE+set}" = set &&
       test -f "${SERPA_PERFORMANCE_FILE}"; then

	# Run with performance check. 
	com="${additional_env} ${program} ${@}"
	echo "${com}"
	${SERPA_GNUTIME} -o $$.time -f "%e %t" ${com} \
	    ${output_redirection_string} ${errors_redirection_string}
	exit_value=${?}

	# Parse temporary file to get recorded time/memory usage
	read etime memory < $$.time
	rm $$.time
	
	# Report the run time/memory usage 
	echo "serpa-perf ${nproc} <${com}> ${etime} ${memory}"

	# Check if time is out of bounds
	search=`grep "serpa-perf ${nproc} <${com}>" ${SERPA_PERFORMANCE_FILE}`
	if test "$?" = "0";then
	    search=`echo $search | sed 's/.*>//'`
	    read compare_time compare_memory <<EOF	    
$search
EOF
        else
	    compare_time=999999
	    compare_memory=0
	fi

	difference=`bc -l <<EOF
a=${etime}
b=${compare_time}
t=${SERPA_DIFFERENCE}
d=a-b
if ( d < 0 ) {
	d = -d
}
p=d/b
print (p < t), "\n"
EOF`
	if test "${difference}" = "1"; then
	    echo "PERFORMANCE PASSED: target ${compare_time} (seconds) current ${etime} sec"
	else
	    echo "PERFORMANCE FAIL: target ${compare_time} (seconds) current ${etime}"
	fi

	# Currently don't do memory check since it is not reporting anything.

    else
	# Run without performance testing
	com="${additional_env} ${program} ${@}"
	echo "${com}"
	eval ${com} ${output_redirection_string} ${errors_redirection_string}
	exit_value=${?}
    fi
  else
    # Because parallel runs sometimes fails due to problems unrelated
    # to the program, we give it several tries before declaring failure.
    if test -z "${SERPA_PARALLEL_TRIES}"; then SERPA_PARALLEL_TRIES=1; fi
    c=1
    while test ${c} -le ${SERPA_PARALLEL_TRIES} ; do
      if test ${c} -gt 1 ; then
        echo "possibly failed.  TRY number ${c}"
      fi
      run_multiproc ${program} ${nproc} "${@}"
      exit_value=${?}
      if [ ${exit_value} = 0 ]; then
        break
      fi
      c=`expr ${c} + 1`
    done
  fi
  # Report pass or fail.
  pf_string=FAILED; test "${exit_value}" = 0 && pf_string=PASSED
  cat <<-_EOM_
	${script_name}::::::::::::::::::::::::::::::::::::::::::::::::
	${pf_string}: ${script_name} ${nproc} ${program} ${@}
	${script_name}================================================
	_EOM_
  # Count failures and possibly exit.
  if test ${exit_value} -ne 0 ; then
    serpa_num_failures=`expr ${serpa_num_failures} + 1`;
    if test ${serpa_num_failures} -eq ${SERPA_MAX_FAILS}; then
      die 1 "FAIL running ${program} with ${nproc} processors"
    fi
  fi
done

if test ${serpa_num_failures} -eq 0; then
  echo "${script_name} ${serpa_args} passed."
else
  echo "${script_name} ${serpa_args} had ${serpa_num_failures} failures."
fi
exit ${serpa_num_failures}

# End of script.
__EOM__
])dnl End of macro to create serpa-run.in


dnl Create the configured file from its .in version.
dnl AC_CONFIG_FILES( btng_serpa_run_fn, chmod +x btng_serpa_run_fn )
dnl IMPORTANT NOTE:
dnl Using the defined symbol btng_serpa_run_fn in AC_CONFIG_FILES
dnl causes automake to complain that btng_serpa_run_fn.in is not found.
dnl Since that symbol is a defined m4 macro, this is an automake deficiency.
dnl Until this is fixed, the above AC_CONFIG_FILES usage will not
dnl work!  You must call AC_CONFIG_FILES with the actual names of
dnl the file instead of the m4-defined symbol btng_serpa_run_fn.



AC_CONFIG_COMMANDS([create-serpa-machine-file],[
# Generate a machine file if needed.
if test ! -r "$serpa_machine_file" ; then
  # Set the serpa machine file name to be in the same directory as
  # the serpa-run script.
  btng_serpa_dir_name=`echo ']btng_serpa_run_fn[' | sed -e ['s:^\([^/]*\)$:./\1:'] -e ['s:/[^/]*$::']`;
  btng_serpa_machine_file="${btng_serpa_dir_name}/serpa.machines"
  hostname="`hostname`"
  cat <<-__EOF__>"$btng_serpa_machine_file"
	$hostname
	$hostname
	$hostname
	$hostname
	__EOF__
fi
],[
dnl Settings to make before running the above.
])


# Add the serpa machine file to the clean list.
btng_serpa_dir_name=`echo ']btng_serpa_run_fn[' | sed -e ['s:^\([^/]*\)$:./\1:'] -e ['s:/[^/]*$::']`;
btng_serpa_machine_file="${btng_serpa_dir_name}/serpa.machines"
DISTCLEANFILES="$DISTCLEANFILES ${btng_serpa_machine_file}"
CASC_AC_LOG_VAR(btng_serpa_machine_file DISTCLEANFILES, in support-serpa-run)


undefine([btng_serpa_run_fn])


# End macro $0
])dnl End CASC_SUPPORT_SERPA macro.
