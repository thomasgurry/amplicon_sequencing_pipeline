from CommLink import *

# Launch plotting listener
inbox_folder = '/home/ubuntu/inbox'
cl = CommLink('plotting')
cl.launch_plotting_listener(inbox_folder)
