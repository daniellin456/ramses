import numpy as np
import struct
import math

# =======================================================================
# Load RAMSES data a la OSIRIS
# =======================================================================
def load_snapshot(nout):
    
    infile = generate_fname(nout)
        
    # Read info file and create info dictionary
    infofile = infile+"/info_"+infile.split("_")[-1]+".txt"
    try:
        with open(infofile) as f:
            content = f.readlines()
        f.close()
    except IOError:
        # Clean exit if the file was not found
        print("Info file not found: "+infofile)
        return 0
    
    info = dict()
    for line in content:
        sp = line.split("=")
        if len(sp) > 1:
            try:
                info[sp[0].strip()] = eval(sp[1].strip())
            except NameError:
                info[sp[0].strip()] = sp[1].strip()
    
    # Read the number of variables from the hydro_file_descriptor.txt
    hydrofile = infile+"/hydro_file_descriptor.txt"
    try:
        with open(hydrofile) as f:
            content = f.readlines()
        f.close()
    except IOError:
        # Clean exit if the file was not found
        print("hydro_file_descriptor not found: "+hydrofile)
        return 0
    # Read the total number of hydro variables
    for line in content:
        sp = line.split("=")
        if len(sp) > 1:
            if sp[0].strip() == "nvar":
                info["nvar"] = int(sp[1].strip())
                break
    # Now go through all the variables and check if they are to be read or skipped
    list_vars = []
    ivar = 0
    for line in content:
        sp = line.split(":")
        if len(sp) > 1:
            list_vars.append(sp[1].strip())
            ivar += 1
            
    # Make sure we always read the coordinates
    list_vars.extend(("level","x","y","z","dx"))
    nvar_read = len(list_vars)
        
    # Now read the amr and hydro files =============================================
    # We have to open the files in binary format, and count all the bytes in the ===
    # file structure to extract just the data we need. =============================
    # See output_amr.f90 and output_hydro.f90 in the RAMSES source. ================
    print("Processing %i files in " % (info["ncpu"]) + infile)
    
    # We will store the cells in a dictionary which we build as we go along.
    # The final concatenation into a single array will be done once at the end.
    data_pieces = dict()
    npieces = 0
    
    # Allocate work arrays
    twotondim = 2**info["ndim"]
    xcent = np.zeros([8,3],dtype=np.float64)
    xg    = np.zeros([info["ngridmax"],3],dtype=np.float64)
    son   = np.zeros([info["ngridmax"],twotondim],dtype=np.int32)
    var   = np.zeros([info["ngridmax"],twotondim,nvar_read],dtype=np.float64)
    xyz   = np.zeros([info["ngridmax"],twotondim,info["ndim"]],dtype=np.float64)
    ref   = np.zeros([info["ngridmax"],twotondim],dtype=np.bool)
    
    iprog = 1
    istep = 10
    ncells_tot = 0
    
    # Loop over the cpus and read the AMR and HYDRO files in binary format
    for k in range(info["ncpu"]):
        
        # Print progress
        percentage = int(float(k)*100.0/float(info["ncpu"]))
        if percentage >= iprog*istep:
            print("%3i%% : read %10i cells" % (percentage,ncells_tot))
            iprog += 1
        
        # Read binary AMR file
        amr_fname = generate_fname(nout,ftype="amr",cpuid=k+1)
        with open(amr_fname, mode='rb') as amr_file: # b is important -> binary
            amrContent = amr_file.read()
        amr_file.close()
        
        # Read binary HYDRO file
        hydro_fname = generate_fname(nout,ftype="hydro",cpuid=k+1)
        with open(hydro_fname, mode='rb') as hydro_file: # b is important -> binary
            hydroContent = hydro_file.read()
        hydro_file.close()
        
        # Need to extract info from the file header on the first loop
        if k == 0:
        
            # nx,ny,nz
            ninteg = 2
            nfloat = 0
            nlines = 2
            nstrin = 0
            nquadr = 0
            offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
            [nx,ny,nz] = struct.unpack("3i", amrContent[offset:offset+12])
            ncoarse = nx*ny*nz
            xbound = [float(int(nx/2)),float(int(ny/2)),float(int(nz/2))]
            
            # nboundary
            ninteg = 7
            nfloat = 0
            nlines = 5
            nstrin = 0
            nquadr = 0
            offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
            nboundary = struct.unpack("i", amrContent[offset:offset+4])[0]
            ngridlevel = np.zeros([info["ncpu"]+nboundary,info["levelmax"]],dtype=np.int32)
            
            # noutput
            ninteg = 9
            nfloat = 1
            nlines = 8
            nstrin = 0
            nquadr = 0
            offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
            noutput = struct.unpack("i", amrContent[offset:offset+4])[0]

        # Read the number of grids
        ninteg = 14+(2*info["ncpu"]*info["levelmax"])
        nfloat = 18+(2*noutput)+(2*info["levelmax"])
        nlines = 21
        nstrin = 0
        nquadr = 0
        offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
        ngridlevel[:info["ncpu"],:] = np.asarray(struct.unpack("%ii"%(info["ncpu"]*info["levelmax"]), amrContent[offset:offset+4*info["ncpu"]*info["levelmax"]])).reshape(info["levelmax"],info["ncpu"]).T
        
        # Read boundary grids if any
        if nboundary > 0:
            ninteg = 14+(3*info["ncpu"]*info["levelmax"])+(10*info["levelmax"])+(2*nboundary*info["levelmax"])
            nfloat = 18+(2*noutput)+(2*info["levelmax"])
            nlines = 25
            nstrin = 0
            nquadr = 0
            offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16 + 4
            ngridlevel[info["ncpu"]:info["ncpu"]+nboundary,:] = np.asarray(struct.unpack("%ii"%(nboundary*info["levelmax"]), amrContent[offset:offset+4*nboundary*info["levelmax"]])).reshape(info["levelmax"],nboundary).T

        # Determine bound key precision
        ninteg = 14+(3*info["ncpu"]*info["levelmax"])+(10*info["levelmax"])+(3*nboundary*info["levelmax"])+5
        nfloat = 18+(2*noutput)+(2*info["levelmax"])
        nlines = 21+2+3*min(1,nboundary)+1+1
        nstrin = 128
        nquadr = 0
        offset = 4*ninteg + 8*(nlines+nfloat) + nstrin + nquadr*16
        key_size = struct.unpack("i", amrContent[offset:offset+4])[0]
        
        # Offset for AMR
        ninteg1 = 14+(3*info["ncpu"]*info["levelmax"])+(10*info["levelmax"])+(3*nboundary*info["levelmax"])+5+3*ncoarse
        nfloat1 = 18+(2*noutput)+(2*info["levelmax"])
        nlines1 = 21+2+3*min(1,nboundary)+1+1+1+3
        nstrin1 = 128 + key_size
        
        # Offset for HYDRO
        ninteg2 = 5
        nfloat2 = 1
        nlines2 = 6
        nstrin2 = 0
            
        # Loop over levels
        for ilevel in range(info["levelmax"]):
            
            # Geometry
            dxcell=0.5**(ilevel+1)
            dx2=0.5*dxcell
            for ind in range(twotondim):
                iz=int((ind)/4)
                iy=int((ind-4*iz)/2)
                ix=int((ind-2*iy-4*iz))
                xcent[ind,0]=(float(ix)-0.5)*dxcell
                xcent[ind,1]=(float(iy)-0.5)*dxcell
                xcent[ind,2]=(float(iz)-0.5)*dxcell
            
            # Cumulative offsets in AMR file
            ninteg_amr = ninteg1
            nfloat_amr = nfloat1
            nlines_amr = nlines1
            nstrin_amr = nstrin1
            
            # Cumulative offsets in HYDRO file
            ninteg_hydro = ninteg2
            nfloat_hydro = nfloat2
            nlines_hydro = nlines2
            nstrin_hydro = nstrin2
                            
            # Loop over domains
            for j in range(nboundary+info["ncpu"]):
                
                ncache = ngridlevel[j,ilevel]
                
                # Skip two lines of integers
                nlines_hydro += 2
                ninteg_hydro += 2
                
                if ncache > 0:
                
                    if j == k:
                        # xg: grid coordinates
                        ninteg = ninteg_amr + ncache*3
                        nfloat = nfloat_amr
                        nlines = nlines_amr + 3
                        nstrin = nstrin_amr
                        for n in range(info["ndim"]):
                            offset = 4*ninteg + 8*(nlines+nfloat+n*(ncache+1)) + nstrin + 4
                            xg[:ncache,n] = struct.unpack("%id"%(ncache), amrContent[offset:offset+8*ncache])
                            
                        # son indices
                        ninteg = ninteg_amr + ncache*(4+2*info["ndim"])
                        nfloat = nfloat_amr + ncache*info["ndim"]
                        nlines = nlines_amr + 4 + 3*info["ndim"]
                        nstrin = nstrin_amr
                        for ind in range(twotondim):
                            offset = 4*(ninteg+ind*ncache) + 8*(nlines+nfloat+ind) + nstrin + 4
                            son[:ncache,ind] = struct.unpack("%ii"%(ncache), amrContent[offset:offset+4*ncache])
                            # var: hydro variables
                            #jvar = 0
                            for ivar in range(info["nvar"]):
                                #if var_read[ivar]:
                                offset = 4*ninteg_hydro + 8*(nlines_hydro+nfloat_hydro+(ind*info["nvar"]+ivar)*(ncache+1)) + nstrin_hydro + 4
                                var[:ncache,ind,ivar] = struct.unpack("%id"%(ncache), hydroContent[offset:offset+8*ncache])
                                #jvar += 1
                            var[:ncache,ind,-5] = float(ilevel+1)
                            for n in range(info["ndim"]):
                                xyz[:ncache,ind,n] = xg[:ncache,n] + xcent[ind,n]-xbound[n]
                                var[:ncache,ind,-4+n] = xyz[:ncache,ind,n]*info["boxlen"]
                            var[:ncache,ind,-1] = dxcell*info["boxlen"]
                            # ref: True if the cell is unrefined
                            ref[:ncache,ind] = np.logical_not(np.logical_and(son[:ncache,ind] > 0,ilevel < info["levelmax"]-1))

                        cube = np.where(ref[:ncache,:])
                        cells = var[cube]
                        ncells = np.shape(cells)[0]
                        if ncells > 0:
                            ncells_tot += ncells
                            npieces += 1
                            # Add the cells in the master dictionary
                            data_pieces["piece"+str(npieces)] = cells
                            
                    # Now increment the offsets while looping through the domains
                    ninteg_amr += ncache*(4+3*twotondim+2*info["ndim"])
                    nfloat_amr += ncache*info["ndim"]
                    nlines_amr += 4 + 3*twotondim + 3*info["ndim"]
                    
                    nfloat_hydro += ncache*twotondim*info["nvar"]
                    nlines_hydro += twotondim*info["nvar"]
            
            # Now increment the offsets while looping through the levels
            ninteg1 = ninteg_amr
            nfloat1 = nfloat_amr
            nlines1 = nlines_amr
            nstrin1 = nstrin_amr
            
            ninteg2 = ninteg_hydro
            nfloat2 = nfloat_hydro
            nlines2 = nlines_hydro
            nstrin2 = nstrin_hydro
    
    # Merge all the data pieces into the master data array
    master_data_array = np.concatenate(list(data_pieces.values()), axis=0)
    
    # Free memory
    del data_pieces,xcent,xg,son,var,xyz,ref
    
    print("Total number of cells loaded: %i" % ncells_tot)
    
    # This is the master data dictionary.
    data = dict()
    for i in range(len(list_vars)):
        theKey = list_vars[i]
        data[theKey] = master_data_array[:,i]

    # Append useful variables to dictionary
    data["unit_d"] = info["unit_d"]
    data["unit_l"] = info["unit_l"]
    data["unit_t"] = info["unit_t"]
    data["boxlen"] = info["boxlen"]
    data["ncells"] = ncells_tot
    data["time"  ] = info["time"]
    data["ngrp"  ] = info["ngrp"]
    
    # Read sink particles if present
    sinkfile = infile+"/sink_"+infile.split("_")[-1]+".csv"
    try:
        sinklist = np.loadtxt(sinkfile,delimiter=",")
        if np.shape(sinklist)[0] == 0:
            data["nsinks"] = 0
        else:
            list_shape = np.shape(np.shape(sinklist))[0]
            if list_shape == 1:
                sinklist = np.reshape(sinklist, (1, np.shape(sinklist)[0]))
                data["nsinks"] = 1
            else:
                data["nsinks"] = np.shape(sinklist)[0]
            try:
                data["r_sink"] = info["ir_cloud"]/(2.0**info["levelmax"])
            except KeyError:
                try:
                    data["r_sink"] = info["ncell_racc"]/(2.0**info["levelmax"])
                except KeyError:
                    data["r_sink"] = 4.0/(2.0**info["levelmax"])
            for i in range(data["nsinks"]):
                data["sink"+str(i+1)] = sinklist[i,:]
    except IOError:
        pass
    
    return data

# =======================================================================
# Generate filename
# =======================================================================
def generate_fname(nout,ftype="",cpuid=1):
        
    number = str(nout).zfill(5)
    infile = "output_"+number
    if len(ftype) > 0:
        infile = infile+"/"+ftype+"_"+number+".out"+str(cpuid).zfill(5)
        
    return infile

# =======================================================================
# Check results against reference solution
# =======================================================================
# - tolerance   : allowed relative difference between the sum over all cells and reference value. Ex: tolerance={"density":1.0e-14}
# - threshold   : relative value below which a vector component is set to zero
# - norm_min    : minimum value for norm, to protect against null vectors
# - min_variance: if the data differs by less than this value from the average value, it is set to the average
def check_solution(data,test_name,tolerance=None,threshold=2.0e-14,norm_min=1.0e-30,min_variance=1.0e-14):

    var_tol = {"all":3.0e-13}
    try:
        for key in tolerance.keys():
            var_tol[key] = tolerance[key]
    except AttributeError:
        pass
    
    # Write dummy file to avoid latex errors
    tex_file = open(test_name+".tex", "w")
    tex_file.write(" \n")
    tex_file.close()
    
    # Find vectors and normalize components
    norms = dict()
    permutations = {"_x":["_y","_z"],"_y":["_x","_z"],"_z":["_x","_y"]}
    for key in sorted(data.keys()):
        norms[key] = 1.0e-10
        if key.endswith("_x") or key.endswith("_y") or key.endswith("_z"):
            rawkey = key[:-2]
            suffix = key[-2:]
            ok = True
            try:
                test = len(data[rawkey+permutations[suffix][0]])
            except KeyError:
                ok = False
            try:
                test = len(data[rawkey+permutations[suffix][1]])
            except KeyError:
                ok = False
            if ok:
                norms[key] = np.sqrt(data[key]**2 + data[rawkey+permutations[suffix][0]]**2 + data[rawkey+permutations[suffix][1]]**2)
                indices = norms[key] < norm_min
                norms[key][indices] = norm_min

    # Compute solution sums
    nvar = len(data.keys())
    sol  = dict()
    ivar = 0
    for key in sorted(data.keys()):
        # Filter out values that are close to the average
        # This is useful if there is a constant non-zero pressure
        # with noise around in the average.
        keyAv = np.average(data[key])
        if keyAv == 0.0:
            keyData = data[key]
        else:
            keyData = np.where(np.abs(data[key]-keyAv)/abs(keyAv) < min_variance,keyAv,data[key])
        # Perform sum
        if key == "density" or \
           key == "thermal_pressure" or \
           key == "total_energy" or \
           key == "temperature" or \
           key.startswith("radiative_energy"):
            solution = np.log10(np.abs(keyData))
        else:
            solution = np.where(np.abs(keyData)<threshold*norms[key],0.0,np.abs(keyData))

        try:
            sol[key] = math.fsum(solution)
        except TypeError:
            sol[key] = solution
    
    ## Overwrite reference solution =====================
    #print("WARNING! Over-writing reference solution")
    #ref_file = open(test_name+"-ref.dat", "w")
    #for key in sorted(data.keys()):
    #    ref_file.write("%s : %.16e\n" % (key,sol[key]))
    #ref_file.close()
    ## ==================================================
    
    # Read reference solution
    ref = dict()
    with open(test_name+"-ref.dat") as f:
        content = f.readlines()
    f.close()
    for line in content:
        sp = line.split(":")
        if len(sp) > 1:
            ref[sp[0].strip()] = eval(sp[1].strip())

    # Checking for errors
    if ref.keys() != sol.keys():
        print("The current and reference solutions do not have the same variables")
        exit()

    ok = True
    # Write error table to tex file
    tex_file = open(test_name+".tex", "w")
    #tex_file.write("\documentclass[12pt]{article}\n")
    #tex_file.write("\usepackage{graphicx,color}\n")
    #tex_file.write("\usepackage[colorlinks=true,linkcolor=blue]{hyperref}\n")
    #tex_file.write("\\begin{document}\n")
    tex_file.write("\\begin{table}[ht]\n")
    tex_file.write("\\scriptsize\n")
    tex_file.write("\\centering\n")
    tex_file.write("\\caption{"+test_name+" error summary}\n")
    tex_file.write("\\begin{tabular}{|l|l|l|l|}\n")
    tex_file.write("\\hline\n")
    tex_file.write("Variable & This run & Reference & Error \\\\\n")
    tex_file.write("\\hline\n")

    # Compute errors
    for key in sorted(ref.keys()):
    #for i in range(np.shape(ref)[0]):
        if sol[key] == ref[key] == 0.0:
            error = 0.0
        else:
            error = abs(sol[key]-ref[key])/min(abs(sol[key]),abs(ref[key]))

        try:
            tol = var_tol[key]
        except KeyError:
            tol = var_tol["all"]

        if error > tol:
            ok = False
            tex_file.write("\\textcolor{red}{%s} & \\textcolor{red}{%.16e} & \\textcolor{red}{%.16e} & \\textcolor{red}{%.16e}\\\\\n" %(key.replace("_"," "),sol[key],ref[key],error))
        else:
            tex_file.write("%s & %.16e & %.16e & %.16e\\\\\n" %(key.replace("_"," "),sol[key],ref[key],error))
    tex_file.write("\\hline\n")
    tex_file.write("\\end{tabular}\n")
    tex_file.write("\\end{table}\n")
    #tex_file.write("\\end{document}\n")
    tex_file.close()

    # Print message if successful
    if ok:
        print("PASSED")
        
    return
