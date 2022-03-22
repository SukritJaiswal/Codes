import pandas as pd

df = pd.read_excel (r'/home/sukrit/Desktop/Project 2020/Spreadsheets/eos_tidal.ods', engine="odf") #place "r" before the path string to address special character, such as '\'. Don't forget to put the file name at the end of the path + '.xlsx'
print (df)

