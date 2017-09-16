Nlmeans, non-local means denoising for Nuke

Plugin by Lewis Saunders (lewis@lewissaunders.com)
Algorithm by: 	Antoni Buades (toni.buades@uib.es
		Bartomeu Coll (tomeu.coll@uib.es)
		Jean-Michel Morel (morel@cmla.ens-cachan.fr)

For more about the algorithm, see:
	http://www.ipol.im/pub/algo/bcm_non_local_means_denoising/

	
Installing
==========
1) Copy the plugin file from linux-32 or mac-intel-32 into ~/.nuke/
2) Add the text from "nlmeans_for_menu.py" to ~/.nuke/menu.py
3) Restart Nuke and you should have a new entry in the Filter menu.


Using
=====
Adjust sigma to remove more or less noise. Increasing the bloc value
sometimes gives better results, but gets slower.

You may want to zoom in, turn on viewer ROI, or crop your image before starting, since it's a pretty slow process.

Enjoy!


License
=======
The original algorithm source code is under the GPLv3 license, which means the plugin is too.
