import socket
maya = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
maya.connect(("127.0.0.1", 5055))
maya.send('unloadPlugin raytracerd_x64.mll')

