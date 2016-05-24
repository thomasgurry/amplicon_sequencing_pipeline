from CommLink import *

# Launch listener
inbox_folder = '/home/ubuntu/inbox'
cl = CommLink('picrust')
cl.launch_picrust_listener(inbox_folder)
