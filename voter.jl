using PyPlot
MersenneTwister()

# Parameters

L = 100

neigh = 0.25
cand  = 0.25
me    = 1-neigh-cand
perc  = 0.1

# Functions

function closest(x,v1,v2)
    if abs(x-v1)/abs(x-v2)<1
        v1
    else
        v2
    end
end

# Start

A = rand(L+2,L+2)

A[1,:] = A[L+1,:]
A[L+2,:] = A[2,:]
A[:,1] = A[:,L+1]
A[:,L+2] = A[:,2]

#plt.clf()
#plt.ion()
#plt.imshow(A[2:L+1,2:L+1])
#plt.show()

v1m = Float64[]
v2m = Float64[]

v1t = Float64[]
v2t = Float64[]

w1t = Float64[]
m2t = Float64[]

    
for k=1:perc*L*L
   push!(v1m,A[rand(2:L+1),rand(2:L+1)])
end
    
for k=1:perc*L*L
   push!(v2m,A[rand(2:L+1),rand(2:L+1)])
end
    
v1 = 0.5*rand()
v2 = 0.5*rand()+0.5

for i=2:L+1
   for j=2:L+1
       if closest(A[i,j],v1,v2)==v1
           B[i,j] = 0
       else
           B[i,j] = 1
       end
   end
end

push!(v1t,v1)
push!(v2t,v2)
push!(w1t,sum(B[2:L+1,2:L+1])/L^2)
push!(m2t,mean(A[2:L+1,2:L+1]))

for i=1:200
    for j=1:100000
    # Pick a site and direction
        x = rand(2:(L+1))
        y = rand(2:(L+1))
        r = rand()
        cc = closest(A[x,y],v1,v2)
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
    v1m = Float64[]
    v2m = Float64[]
    
    for k=1:perc*L*L
       push!(v1m,A[rand(2:L+1),rand(2:L+1)])
    end
    
    for k=1:perc*L*L
       push!(v2m,A[rand(2:L+1),rand(2:L+1)])
    end

    v1 = (v1+median(v1m))/2
    v2 = (v2+median(v2m))/2

    for i=2:L+1
    for j=2:L+1
        if closest(A[i,j],v1,v2)==v1
            B[i,j] = 0
        else
            B[i,j] = 1
        end
    end
    end

    push!(v1t,v1)
    push!(v2t,v2)
    push!(w1t,sum(B[2:L+1,2:L+1])/L^2)
    push!(m2t,mean(A[2:L+1,2:L+1]))
end

plt.close()
plot(v1t)
plot(v2t)
plot(w1t)
plot(m2t)

[v1t v2t w1t m2t]


