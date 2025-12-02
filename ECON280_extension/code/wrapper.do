/*=======================================

* Author: Zaria Roller
* Date: 12/5/25

* Extension of "Childcare, Labor Supply, and 
Business Development: Experimental Evidence 
from Uganda" 

* This file runs all the code for my extension. 

========================================*/

clear all

* set globals

global dir = "/Users/zariaroller/Downloads/ECON280_extension"

global code  		"$dir/code"
global data  		"$dir/data/raw"
global constructed 	"$dir/data/constructed"
global results 		"$dir/results"

* run analysis file
do analysis.do

* create tables
do tables.do
