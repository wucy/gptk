#
# GPTK
# Xerces C library installer
#
# (c) Remi Barillec, 2009 (installer)
#

XERCES_URL=http://mirror.fubra.com/ftp.apache.org/xerces/c/3/binaries
XERCES_ARCHIVE=xerces-c-3.0.1-x86-linux-gcc-3.4.tar.gz
ROOT=.

# Install xerces
echo "*********************************************************************"
echo "Installing XERCES C library"
echo "*********************************************************************"
echo "Downloading archive..."
wget -nv $XERCES_URL/$XERCES_ARCHIVE
#
echo "Extracting library..."
tar -xzf $XERCES_ARCHIVE -C $ROOT/lib
#
echo "Deleting archive..."
rm $XERCES_ARCHIVE
#
echo "Creating symlink..."
cd $ROOT/lib
ln -s xerces-c-* xercesc
cd ..
#
echo "Done"


