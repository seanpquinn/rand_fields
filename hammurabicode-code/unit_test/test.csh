#!/bin/csh


echo ""
echo "-----------"
echo "Running hammurabi unit test"
echo ""

if (! -d out) mkdir out
rm out/*

echo "Running hammurabi for a map"
../run/hammurabi inputs/params.txt >& out/test.log 

if ($? == 0) then 
echo "  Hammurabi ran successfully." 
else 
echo "  Hammurabi has hit an error:"
tail out/test.log
endif

echo "Running hammurabi with a random component and Faraday rotation at highest res, four shells, on a patch of sky"
../run/hammurabi inputs/params8.txt >& out/test8.log

if ($? == 0) then
echo "  Hammurabi ran successfully."
else
echo "  Hammurabi has hit an error:"
tail out/test8.log
endif

echo ""
echo "Diffing lists (you should see nothing):"
diff ref/test8_list.txt out/test8_list.txt | tail
echo ""
echo "Done"

echo "**NOW COMPARE map unit_test/ref/test_obs.fits with unit_test/out/test_obs.fits"
echo ""
