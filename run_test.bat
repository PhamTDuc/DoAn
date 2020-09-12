@echo off 
echo ---------------Start Running Test for DoAn---------------------
@python -m unittest discover -v 4 -s tests -t . 
@pause