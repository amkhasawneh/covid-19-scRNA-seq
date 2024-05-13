#!/usr/bin/bash


try=0
COMPLETE_CONDITION=0

echo "START"

until [ "$lastresult" = "$COMPLETE_CONDITION" ]; do
  let "try+=1"
  echo "Try $try ..."
  ncftpput -B 33554432 -z -u 'geoftp' -p 'inAlwokhodAbnib5' -v -R ftp-private.ncbi.nlm.nih.gov ./uploads/amkhasawneh_iTgdkk8y GEO_COVID-19


  let "lastresult=$?"
  echo "Last Resultcode: $lastresult"
done

echo "UPLOAD COMPLETED AFTER $try TRY(S)"

exit 0
