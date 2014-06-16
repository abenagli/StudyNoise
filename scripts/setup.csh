if (${?STUDYNOISE}) then
echo "already set"

else
setenv THISDIR `pwd`
setenv STUDYNOISE ${THISDIR}/
setenv PATH ${PATH}:${THISDIR}/bin
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${THISDIR}/lib
endif
