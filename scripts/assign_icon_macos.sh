#!/bin/sh
# Create temporary resource file which points to tmpicons.icns:
echo "read 'icns' (-16455) \"assets/icons/chaos.icns\";" >> chaos.rsrc

# append this resource to the file you want to icon-ize.
Rez -a chaos.rsrc -o bin/main

# Use the resource to set the icon.
SetFile -a C bin/main

# Clean up
rm chaos.rsrc