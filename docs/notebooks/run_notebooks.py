# https://nbconvert.readthedocs.io/en/latest/usage.html#converting-multiple-notebooks
c = get_config()
c.NbConvertApp.notebooks = [
	"docs/notebooks/Visualization.ipynb",
	"docs/notebooks/Kamodo.ipynb",
	"docs/notebooks/Kamodofying_Models.ipynb",
	]