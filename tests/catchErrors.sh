fail=$(echo "${1}" | grep -E 'FAILED: [1-9]' | wc -l)
fail2=$(echo "${1}" | grep -E 'Error: Command failed' | wc -l)
err=$(echo "${1}" | grep -E '[1-9] error' | wc -l)
err2=$(echo "${1}" | grep -E '[1-9] ERROR' | wc -l)
if [ ${fail} -ne 0  ] || [ ${fail2} -ne 0 ]; then
	echo "Failure occurred!"
	exit 229
elif [ ${err} -ne 0 ] || [ ${err2} -ne 0 ]; then
	echo "Error occurred!"
	exit 230
else
	echo "No errors or failures in this test!"
fi
