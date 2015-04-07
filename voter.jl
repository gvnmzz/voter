using PyPlot
MersenneTwister()

# Parameters
nc = 3                  #number of candidates
L = 100                 #lattice dimension
nswps = 200             #number of sweeps

neigh = 0.5             #neighbor weight
cand  = 0.0             #closest candidate weight
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

# Start

A = rand(L+2,L+2)
B = copy(A)

A[1,:] = A[L+1,:]
A[L+2,:] = A[2,:]
A[:,1] = A[:,L+1]
A[:,L+2] = A[:,2]

plt.clf()
plt.ion()
plt.imshow(A[2:L+1,2:L+1])
plt.show()

vm = zeros(Float64,int(perc*L*L),nc)
vt = zeros(Float64,nswps+1,nc)
wt = zeros(Float64,nswps+1)


    
for k=1:perc*L*L
  for l=1:nc
   vm[k,l]=A[rand(2:L+1),rand(2:L+1)]
  end
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
    
    plt.imshow(A[2:L+1,2:L+1])
    plt.draw()
    
# Change candidates
    vm = zeros(Float64,int(perc*L*L),nc)
    
    for k=1:perc*L*L
    for l=1:nc
      vm[k,l]=A[rand(2:L+1),rand(2:L+1)]
    end
    end
    
    for k=1:nc
      v[k] = (v[k]+median(vm[k,:]))/2
    end

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

    vt[i,:] =  v'
    wt[i] = sum(B[2:L+1,2:L+1])/L^2
end

outfile = open("voter_mod.txt", "w")
    writedlm(outfile,A[2:L+1,2:L+1],"\t")
close(outfile)


plt.close()
plot(vt[:,1])
plot(vt[:,2])
plot(vt[:,3])
plot(wt)

[vt wt]


