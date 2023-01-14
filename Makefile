website:
	emacs -Q --script publish.el

publish: website
	rsync -e ssh -vr docs/ ghislain@fdb:/home/www/ecology.ghislainv.fr/

clean:
	rm -r docs/*
