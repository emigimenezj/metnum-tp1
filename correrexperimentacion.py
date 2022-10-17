import os
import sys

if(len(sys.argv)<2):
	print("Usar argumento run para correr la experimentación y clear para limpiar el output de los notebooks.")

elif (sys.argv[1] == "run"):
	os.system("jupyter nbconvert --execute --to notebook --inplace Experimentacion\ parte1.ipynb ")
	os.system("jupyter nbconvert --execute --to notebook --inplace Experimentacion\ parte2.ipynb ")

elif (sys.argv[1] == "clean"):
	os.system("jupyter nbconvert --clear-output --inplace Experimentacion\ parte1.ipynb")
	os.system("jupyter nbconvert --clear-output --inplace Experimentacion\ parte2.ipynb")

else:
	print("Usar argumento run para correr la experimentación y clear para limpiar el output de los notebooks.")
