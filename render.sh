ffmpeg -framerate 144 -i frames/frame_%05d.svg -c:v mpeg4 -q:v 4 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p result.mp4
rm frames/*.svg

