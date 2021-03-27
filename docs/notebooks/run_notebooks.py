# https://nbconvert.readthedocs.io/en/latest/usage.html#converting-multiple-notebooks
c = get_config()
c.NbConvertApp.notebooks = [
	"Visualization.ipynb",
	]