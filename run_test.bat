@echo off
cls
@rem -s for print STDOUT, -v for VERBOSE
echo $E[91m---------------Start Running Test for DoAn---------------------$E[0m
@py -3.8 -m pytest tests/ -v -s
echo ---------------End Test for DoAn-------------------------------
@pause
