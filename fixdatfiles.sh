#!/bin/bash
string='mmol/g'
newstr='loading(mmol/g)'



sed -i "s|$string|$newstr|g" *313K.dat
#sed -i 's/,/ /g' *.csv
