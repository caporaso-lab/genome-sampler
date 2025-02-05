html:
	jupyter book build --html
	# cp -r data/ _build/html/data/

serve:
	npx serve _build/html/ -p 4000

clean:
	rm -rf _build/html/