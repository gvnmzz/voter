using PyPlot
MersenneTwister()

# Parameters
nc = 2                  #number of candidates
L = 100                 #lattice dimension
nswps = 200             #number of sweeps

neigh = 0.3             #neighbor weight
cand  = 0.1             #closest candidate weight
me    = 1-neigh-cand    #own weight
perc  = 0.1             #percentage of data candidates have

# Functions

function closest(x,v)
    res = v[1]
    for i = 2:length(v)
    if abs(x-res) < abs(x-v[i])
        res = v[i]
    end
    end
    res
end
outfile = open("voter_m2.txt", "w")

# Start
for perc2 = 0.01:0.002:0.2
A = rand(L+2,L+2)
B = copy(A)

A[1,:] = A[L+1,:]
A[L+2,:] = A[2,:]
A[:,1] = A[:,L+1]
A[:,L+2] = A[:,2]

#plt.clf()
#plt.ion()
#plt.imshow(A[2:L+1,2:L+1])
#plt.show()

vm = zeros(Float64,int(perc*L*L))
vm2 = zeros(Float64,int(perc2*L*L))
vt = zeros(Float64,nswps,nc)
wt = zeros(Float64,nswps,nc)


    
for k=1:perc*L*L
   vm[k]=A[rand(2:L+1),rand(2:L+1)]
end
for k=1:perc2*L*L
   vm2[k]=A[rand(2:L+1),rand(2:L+1)]
end   
    
v = rand(nc)

for i=2:L+1
   for j=2:L+1
       if closest(A[i,j],v)==v[1]
           B[i,j] = 0
       elseif closest(A[i,j],v)==v[2]
           B[i,j] = 1
       else
           B[i,j] = 2
       end
   end
end

vt[1,:] =  v'
wt[1] = sum(B[2:L+1,2:L+1])/L^2

for i=1:nswps
    for j=1:10000
    # Pick a site and direction
        x = rand(2:(L+1))
        y = rand(2:(L+1))
        r = rand()
        cc = closest(A[x,y],v)
    # Average with your neighbor
        if r<0.25
            A[x,y] = me*A[x-1,y]+neigh*A[x,y]+cand*cc
        elseif r>0.25 && r<0.5
            A[x,y] = me*A[x,y+1]+neigh*A[x,y]+cand*cc
        elseif r>0.5 && r<0.75
            A[x,y] = me*A[x+1,y]+neigh*A[x,y]+cand*cc
        elseif r>0.75
            A[x,y] = me*A[x,y-1]+neigh*A[x,y]+cand*cc
        end
        
    # Update boundaries        
        if x==2
            A[L+2,y]=A[x,y]
        end
        if y==2
            A[x,L+2]=A[x,y]
        end
        if x==L+1
            A[1,y]=A[x,y]
        end
        if y==L+1
            A[x,1]=A[x,y]
        end
    # End update boundaries         
    end

#    plt.imshow(A[2:L+1,2:L+1])
#    plt.draw()

# Change candidates
    vm = zeros(Float64,int(perc*L*L))
    vm2 = zeros(Float64,int(perc2*L*L))
    
    for k=1:perc*L*L
        vm[k]=A[rand(2:L+1),rand(2:L+1)]
    end
    for k=1:perc2*L*L
        vm2[k]=A[rand(2:L+1),rand(2:L+1)]
    end   

      v[1] = (v[1]+median(vm))/2
      v[2] = (v[2]+median(vm2))/2
      
    for l=2:L+1
    for j=2:L+1
       if closest(A[l,j],v)==v[1]
           B[l,j] = 0
           wt[i,1] += 1
       elseif closest(A[l,j],v)==v[2]
           B[l,j] = 1
           wt[i,2] += 1
       else
           B[l,j] = 2
           wt[i,3] += 1
       end
    end
    end
    vt[i,:] =  v'

end


#    writedlm(outfile,[vt wt],"\t")
    writedlm(outfile,[perc2 mean(wt[100:200,1]) std(wt[100:200,1])],"\t")

end
close(outfile)


#plt.close()
#plot(vt[:,1])
#plot(vt[:,2])
#plot(vt[:,3])
#plot(wt/2)

#[vt wt]
#println(mean(wt[50:end]/2)," ",std(wt[50:end])/2)

