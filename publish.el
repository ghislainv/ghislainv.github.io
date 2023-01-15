;;; publish.el --- generate and publish my blog -*- lexical-binding: t -*-
;;; Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
;;; Commentary:
;;; Documentation at https://orgmode.org/worg/org-tutorials/org-publish-html-tutorial.html
;;; Code:

;; Package installation
(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(add-to-list 'package-archives '("gnu" . "https://elpa.gnu.org/packages/") t)
(package-initialize)

;; Fetch the list of packages available
(unless package-archive-contents
  (package-refresh-contents))

;; Install dependencies
(defvar my-package-list)
(setq my-package-list '(htmlize))
(dolist (i-package my-package-list)
  (unless (package-installed-p i-package)
    (package-install i-package)))

;; Load the publishing system
(require 'ox-publish)

(defun my/relative-path-expand (path)
  "Expand relative PATH from current buffer or file to a full path."
  (concat
   (if load-file-name
       (file-name-directory load-file-name)
     default-directory)
   path))

(defvar mywebsite-base-directory
  (my/relative-path-expand "org/")
  "The `base-directory' for mywebsite project.")

(defvar mywebsite-publish-directory
  (my/relative-path-expand "docs/")
  ;; "/ssh:ghislain@fdb:/home/www/mywebsite.org/"
  "The `publishing-directory' for mywebsite project.")

;; Don't show section numbers or table of contents by default
(setq org-export-with-section-numbers nil
      org-export-with-toc             nil)

;; Link abbreviation
;; https://orgmode.org/manual/Link-Abbreviations.html
(setq org-link-abbrev-alist
      '(("doi" . "https://doi.org/")))

;; HTML
(require 'ox-html)

;; Set preambule
(setq org-html-preamble t
      org-html-preamble-format
       `(("en" ,(with-temp-buffer
		  (insert-file-contents-literally "./org/resources/preamble.html")
		  (buffer-substring-no-properties (point-min) (point-max))))))

;; Set postambule
(setq org-html-postamble t
      org-html-postamble-format
       `(("en" ,(with-temp-buffer
		  (insert-file-contents-literally "./org/resources/postamble.html")
		  (buffer-substring-no-properties (point-min) (point-max))))))

;; Enable HTML5
(setq org-html-html5-fancy t
      org-html-doctype     "html5")

;; Disable ox-html's default CSS and JavaScript
(setq org-html-head-include-default-style nil
      org-html-head-include-scripts       nil)

;; Render ~verbatim~ as kbd tag in HTML
(add-to-list 'org-html-text-markup-alist '(verbatim . "<kbd>%s</kbd>"))

;; Enable code syntax highlighting if available
(if (package-installed-p 'htmlize)
    (setq org-html-htmlize-output-type 'css)
  (setq org-html-htmlize-output-type 'nil))

(setq org-publish-project-alist
      `(("pages"
	 :base-directory ,mywebsite-base-directory
	 :base-extension "org"
	 :exclude "setup.org"
	 :publishing-directory ,mywebsite-publish-directory
	 :publishing-function org-html-publish-to-html
	 )
	("blog"
	 :base-directory ,(concat mywebsite-base-directory "blog/")
	 :base-extension "org"
	 :publishing-directory ,(concat mywebsite-publish-directory "blog/")
	 :publishing-function org-html-publish-to-html
	 :auto-sitemap t
	 :sitemap-title "Blog Posts"
	 :sitemap-filename "index.org"
	 :sitemap-sort-files anti-chronologically
	 )
	("static"
	 :base-directory ,mywebsite-base-directory
	 :base-extension "ico\\|html\\|css\\|scss\\|woff2\\|jpg\\|gif\\|png\\|txt\\|pdf\\|zip\\|kmz\\|R\\|docx"
	 :recursive t
	 :include ("CNAME" ".nojekyll" "keybase.txt")
	 :publishing-directory ,mywebsite-publish-directory
	 :publishing-function org-publish-attachment)
	("mywebsite" :components ("pages" "blog" "static"))))

;; Uncomment to force full site regeneration
(org-publish "mywebsite" t)

;; Message
(message "Build complete!")

;;; publish.el ends here
