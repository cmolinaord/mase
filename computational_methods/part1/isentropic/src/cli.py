def input(w, opt):
	# Printing input data
	print("Maximum iterations = %i" % opt.itermax)
	print("L = %1.1f" % w.L)
	print("H = %1.1f" % w.H)
	print("[nx,ny] = [%i,%i]" % (w.nx,w.ny))
	print("Precission = %1.1E" % opt.precission)
	print("Compressible? %s" % ('Yes' if opt.compressible else 'No'))
	print("###########################")
	print(" ")
