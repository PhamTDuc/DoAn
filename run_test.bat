rem -s for print STDOUT, -v for VERBOSE
@echo off 
echo ---------------Start Running Test for DoAn---------------------
@py -3.9 -m pytest tests/ -v -s
echo ---------------End Test for DoAn-------------------------------
@pause
