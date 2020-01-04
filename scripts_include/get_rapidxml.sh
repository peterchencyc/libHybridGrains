#!/bin/bash

actual_rapidxml_zip_md5="7b4b42c9331c90aded23bb55dc725d6a"
rapidxml_url="http://downloads.sourceforge.net/project/rapidxml/rapidxml/rapidxml%201.13/rapidxml-1.13.zip"
rapidxml_file_name="rapidxml-1.13.zip"
extracted_rapidxml_name="rapidxml-1.13"
# md5 on installed RapidXml source files
actual_installed_rapidxml_md5="d698fa3e3547874766de5ba8cb2177c9"

# Verify that curl is installed
command -v curl >/dev/null 2>&1 || { echo >&2 "Error, please install curl and rerun get_eigen.sh."; exit 1; }

# Verify that md5sum is installed
command -v md5sum >/dev/null 2>&1 || command -v md5 >/dev/null 2>&1 || { echo >&2 "Error, please install md5 or md5sum and rerun get_dependencies.sh."; exit 1; }
if command -v md5 > /dev/null 2>&1 ; then
  md5sum()
  {
    md5 "$@" | sed -e 's#^MD5 [(]\(.*\)[)] = \(.*\)$#\2 \1#'
  }
  export -f md5sum
fi

# Verify that unzip is installed
command -v unzip >/dev/null 2>&1 || { echo >&2 "Error, please install unzip and rerun get_rapidxml.sh."; exit 1; }

# If the output directory exists
if [ -d "include/rapidxml" ]; then
  # If the RapidXml install is up to date, no action is needed
  computed_installed_rapidxml_md5=`find include/rapidxml -type f -name '*.hpp' -exec bash -c 'md5sum "$0" "$@"' {} + | awk '{print $2$1}' | sort -fd | md5sum | cut -c -32`
  if [ "$computed_installed_rapidxml_md5" == "$actual_installed_rapidxml_md5" ]
  then
    echo "RapidXml library is already up to date, no further action is needed."
    exit 0
  fi
  # Otherwise, the checksum is incorrect, warn the user and exit
  echo "Error, directory include/rapidxml has an incorrect checksum, please run remove_rapidxml.sh and rerun get_rapidxml.sh."
  exit 1
fi

echo "Installing RapidXml"

# Create a temporary working directory
temp_dir_name=`uuidgen`
if [ -d "$temp_dir_name" ]; then
  echo "Error, temporary working directory $temp_dir_name exists, this is a bug. Please contact the maintainer."
  exit 1
fi
echo "--->  Creating temporary directory $temp_dir_name"
mkdir $temp_dir_name

function cleanup {
  echo "--->  Removing temporary directory $temp_dir_name"
  rm -fr "$temp_dir_name"
}
trap cleanup EXIT

# Download RapidXml
echo "--->  Downloading RapidXml source"
curl -s -L -o "$temp_dir_name/$rapidxml_file_name" "$rapidxml_url"
if [ $? -ne 0 ]
then
  echo "Error, failed to download RapidXml from $rapidxml_url."
  exit 1
fi

# Run a checksum on the download
echo "--->  Verifying RapidXml checksum"
computed_rapidxml_zip_md5=`md5sum $temp_dir_name/$rapidxml_file_name | cut -c -32`
if [ "$actual_rapidxml_zip_md5" != "$computed_rapidxml_zip_md5" ]
then
  echo "Error, md5 checksum for $rapidxml_file_name does not match $actual_rapidxml_zip_md5."
  exit 1
fi

# Extract the zip archive
echo "--->  Extracting RapidXml"
unzip -q "$temp_dir_name"/"$rapidxml_file_name" -d "$temp_dir_name"
# Move the source to its final location
echo "--->  Moving RapidXml to destination"
mkdir -p include/rapidxml
mv $temp_dir_name/$extracted_rapidxml_name/*.hpp include/rapidxml/

trap - EXIT
cleanup
echo "Successfully installed RapidXml"
