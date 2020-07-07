pyinstaller pyTiffMicasense.py \
  --add-binary "$(command -v exiftool)":"." \
  --hidden-import='pkg_resources.py2_warn' \
  --hidden-import='pyproj.datadir' \
  --onefile \
  --clean \
  --windowed

#pyinstaller pyTiff.py \
#  --add-binary "$(command -v exiftool)":"." \
#  --hidden-import='pkg_resources.py2_warn' \
#  --hidden-import='pyproj.datadir' \
#  --onefile \
#  --clean