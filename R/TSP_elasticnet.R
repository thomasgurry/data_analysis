
#Generate N cities randomly in a unit square. 
random.cities <- function(N){ 

	points = matrix(runif(2*N,min=0,max=1),nrow=N,ncol=2) 

} 


#Function to calculate Euclidean distance between two points. 
#’y’ may be a matrix of points. 

distance <- function(x,y){ 

	if(length(y)>2){ 
		dimension = dim(y)[1] 
	}else{ 
		dimension = 1 
	} 

	dist = rep(0,dimension) 

	if(dimension>1){ 
		for(i in 1:dimension){ 
			dist[i] = sqrt(sum((x-y[i,])^2)) 
		} 
	}else{ 
		for(i in 1:dimension){ 
			dist[i] = sqrt(sum((x-y[i])^2)) 
		} 
	} 

return(dist) 

} 


#Function phi(d,K). 
phi.function <- function(d,K){ 
	res = exp(-(d^2)/(2*(K^2))) 
	return(res) 
}	 

#Implementation of elastic net for a set of N points. 

TSP <- function(points,N){ 
	#par(ask=TRUE) 
	centroid = c(sum(points[,1])/N,sum(points[,2])/N) 
				
	M = 2.5*N 
	path.points = matrix(0,nrow=M,ncol=2)
	alpha = 0.2 
	beta = 2.0 
	K = 0.2 
	weights = matrix(0,nrow=N,ncol=M) 
	delta.y = matrix(1,nrow=M,ncol=2) 
	x.max = max(max(points[,1]),max(path.points[,1])) 
	x.min = min(min(points[,1]),min(path.points[,1])) 
	y.max = max(max(points[,2]),max(path.points[,2])) 
	y.min = min(min(points[,2]),min(path.points[,2])) 
	range = mean((x.max-x.min),(y.max-y.min)) 
	print(range) 
	
	#Set M path points randomly on a circle of radius range/6. 
	initial.radius = range/10
	which.angles = seq((2*pi)/M,2*pi,by=(2*pi)/M) 
	path.points[,1] = cos(which.angles)*initial.radius + centroid[1] 
	path.points[,2] = sin(which.angles)*initial.radius + centroid[2] 

	iteration = 0 
	change = 1 
	plot.path(points,path.points,centroid,K) 
	while(iteration<5000){ 
		if(K<(0.05+0.0001) & K>(0.05-0.0001)){ 
			plot.path(points,path.points,centroid,K) 
		} 
		if(K<(0.01+0.0001) & K>(0.01-0.0001)){ 
			plot.path(points,path.points,centroid,K) 
		} 

		if(iteration==4999) 
			plot.path(points,path.points,centroid,K) 

		#Reduce value of K by 1% every 5 iterations. 
		iteration = iteration + 1 
		#print(iteration) 

		if(iteration%%5 == 0) 
			K = 0.99*K 
			
		#Compute w_ij according to (2) in Durbin & Willshaw. 
		for(j in 1:M){ 
			for(i in 1:N){ 
				denominator = 0 
				dist1 = sqrt((points[i,1]-path.points[,1])^2 + (points[i,2]-path.points[,2])^2) 
				denominator = denominator + sum(exp(-(dist1^2)/(2*K^2))) 
				if(denominator==0) 
					denominator= 10^(-10) 
				dist = sqrt((points[i,1]-path.points[j,1])^2 + (points[i,2]-path.points[j,2])^2) 
				weights[i,j] = exp(-(dist^2)/(2*K^2)) 
				weights[i,j] = weights[i,j]/denominator 
			} 
		} 

		#Compute delta y_j. 
		for(j in 1:M){ 
			#Include periodic boundary conditions, for a = j+1, j, k = j-1 
			if(j==1){ 
				k = M 
			}else{ 
				k = j-1 
			} 
			if(j==M){ 
				a = 1 
			}else{ 
				a = j+1 
			} 
			
			term = 0 
			for(n in 1:N){ 
				term = term + weights[n,j]*(points[n,]-path.points[j,]) 
			} 

			delta.y[j,1] = alpha*term[1] + beta*K*(path.points[a,1]-2*path.points[j,1]+path.points[k,1]) 
			delta.y[j,2] = alpha*term[2] + beta*K*(path.points[a,2]-2*path.points[j,2]+path.points[k,2]) 
		} 
		
		path.points = path.points + delta.y 
		plot.path(points,path.points,centroid,K) 
	} 

	path.length = 0 
	for(x in 1:M){ 
		if(x==M){ 
			y = 1 
		}else{ 
			y = x+1 
		} 
		path.length = path.length + sqrt(sum(path.points[x,]-path.points[y,])^2) 
	} 

print(paste("Path length: ",path.length)) 
return(path.points) 
} 


#Plot points. 
plot.path <- function(points,path.points,centroid,K){ 
	plot(points[,1],points[,2],type="p",col="black",main=paste("K = ",K)) 
	points(path.points,col="red") 
	centroid = rbind(centroid,centroid) 
	points(centroid,col="blue") 
	lines(path.points,col="red") 
	lines(rbind(path.points[1,],path.points[dim(path.points)[1],]),col="red") 
} 

#Generate points and run TSP. 
go.TSP <- function(N){ 
	if(burma==1){ 
		points = cities[,2:3] 
	}else{ 
		points = random.cities(N) 
	} 
	print(points) 
	path.points = TSP(points,N) 
} 


#######################################
#  CARE: RUN THIS WITH mod(N,10)==0   #    
#######################################

letsgo <- function(N){

	points = random.cities(N)
	path.points = TSP(points,N)

}
