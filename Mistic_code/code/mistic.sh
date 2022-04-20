
find . | grep .git | xargs rm -rf
find . -name ".DS_Store" -delete
bokeh serve --port 5098 --show image_tSNE_GUI
