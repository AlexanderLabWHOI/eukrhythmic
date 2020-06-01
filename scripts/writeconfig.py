## SIMPLE SCRIPT FOR ALLOWING USER TO SPECIFY CONFIG ENTRIES INTERACTIVELY
import os
import sys
import yaml

print('Would you like to input setup values interactively (y/n)?')
x = input()
if (x != "y") & (x != "yes"):
    sys.exit(0)

with open('config.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
print("Simply click enter if you do not have input for any of the statements.")
print('What directory are your input files located in?')
x = input()