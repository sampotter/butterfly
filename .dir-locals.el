;;; Directory Local Variables            -*- no-byte-compile: t -*-
;;; For more information see (info "(emacs) Directory Variables")

((nil . ((eval . (define-skeleton
				   bf-error-skeleton
				   "Inserts an error handling section delimited by
 calls to BF_ERROR_BEGIN and BF_ERROR_END."
				   nil
				   \n > "BF_ERROR_BEGIN();"
				   \n >
				   \n > _
				   \n >
				   \n > "BF_ERROR_END() {"
				   \n > "BF_DIE();"
				   \n > "}" >))
		 (projectile-indexing-method . native)))
 (c++-mode . ((indent-tabs-mode . nil)
			  (c-basic-offset . 2)
			  (comment-style . multi-line)))
 (c-mode . ((indent-tabs-mode . nil)
			(c-basic-offset . 2)
			(comment-style . multi-line)
			(visual-line-mode . t)
			(visual-fill-column-mode . nil)
			(adaptive-wrap-prefix-mode . t)
			(adaptive-wrap-extra-indent . 4)))
 (python-mode . ((subword-mode)))
 (cython-mode . ((subword-mode))))
