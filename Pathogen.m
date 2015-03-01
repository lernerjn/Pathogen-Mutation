function [tout, numberOfInfected, PopWhoToInfect, time] = Pathogen(N, tauR, tauI, K, hthr, virusSize, mu ,p, infectProb, tend) 
% Author: Jeremy Lerner, University of Arizona 
%
% This program models the spread of a pathogen that mutates around a population: 
%
% Outputs: [tout, numberOfInfected, PopWhoToInfect, time]
%       tout: simply an array with 1:1:tend, to plot against the number of
%       infected
%
%       numberOfInfected: the number of people infected at each time step,
%       so numberOfInfected(t) is the number of sick people at time t
%
%       PopWhoToInfect: a representation of the small world network, that
%       is PopWhoToInfect(i,:), for a person with no random connections
%       would be [i - K , i - K + 1, ... , i + 1 , ..., i + K], and with
%       increasing p, there is an inreasing chance of some of the entries
%       will be a random connection
%
%       time: the number of seconds the program took to run
%       
% Inputs: Pathogen(N, tauR, tauI, K, hthr, virusSize mu ,p, infectProb, tend)
%                  1    2     3   4    5        6    7   8     9         10
% 
%       N:  size of the population
%
%       tauR: immunity duration (for example, if you  get a virus at 
%       time 0, then at time tauR you are no longer immune to that
%       virus or viruses similar to that one)
%
%       tauI: infection duration
%
%       K: K is defined such that 2 * K is the number of connections per person
%
%       hthr is the hamming distance threshold, so if a virus is different than
%       the viruses in a person's history by more than the hamming distance, then
%       the person is not immune
%       Example: if the virus is [0 0 0 0] and a person has had virus [0 1 1 1],
%       then the hamming distance is 3, because 3 bits are different.
%
%       virusSize: the number of bits in the virus bit string
%
%       mu: the probability that the virus mutates at each transmission
%
%       p: the probability of a random connection in the network, each
%       connection has a p chance of being replaced with a random
%       connection
%
%       infectProb: the probabilty that the virus infects a person if they come
%       in contact with an infected person
%
%       tend: the number of time steps the model will run for.
%
% An interesting case, these conditions cause a damped oscillation
% [toutg, numberOfInfectedg, PopWhoToInfect]= Pathogen(1e5, 150, 0, 2, 2, 10, 0.01 ,0.01, 1, 2000);
tic;

VH = struct([]); %setting up the virus history structure array

for n = 1 : N
    VH(n).viruses = []; %these will be filled with the viruses that people have had
    VH(n).times = []; %these are the times they go those viruses
    VH(n).virus = []; %this is the virus the person is currently infected with, if any.
    VH(n).time = []; %How long the person has been infected with the current virus they have
end


Pop = zeros(N,1); % 0 means the person is susceptible and 1 means they are infected

%Infecting 1 random person at the beginning
personInfected = unidrnd(N);
Pop( personInfected ) = 1;
VH(personInfected).virus = zeros(1,virusSize);
VH(personInfected).times  = 0;
VH(personInfected).viruses = zeros(1,virusSize);
VH(personInfected).time = 0;


PopWhoToInfect = zeros(N , 2 * K );

%Building the network:

%setting up the edge cases on the network with wraparound connections

%the right side edge case
for i = N - K + 1 : N 
    PopWhoToInfect(i,1:K) = (i-K):(i-1);
    PopWhoToInfect(i,K+1:2*K) = [(i+1:N) 1:(K-(N-i))];
    
    for ( connectR = 1: 2 * K)
    
        if (rand(1) < p )
            person  = unidrnd(N);
            while ( ~(isempty(find(PopWhoToInfect(i,:) == person))) || person == i )
                person = unidrnd(N);
            end

            PopWhoToInfect(i, connectR  ) = unidrnd(N);


        end
    end
end 

%the left side edge case
for i = 1 : K
    PopWhoToInfect(i,K+1:2*K) =  i+1:(K+i);
    PopWhoToInfect(i,1:K) =  [(N - (K - i)):N  1:i-1];
    
    for (connectL = 1: 2*K)
    
        if (rand(1) < p )
            person  = unidrnd(N);
            while ( ~(isempty(find(PopWhoToInfect(i,:) == person))) || person == i )
                person = unidrnd(N);
            end

            PopWhoToInfect(i, connectL ) = unidrnd(N);


        end
    end

end

    
%the main part of the network
for index = K + 1 : N  -  K 
    
    PopWhoToInfect(index,1:K) = (index-K):(index-1);
    PopWhoToInfect(index,K+1:2*K) =  index+1:(K+index);
    
    PopWhoToInfect(index, 2 * K ) = index + K;

    
    for (connectM = 1:2*K)

        if (rand(1) < p ) %the random connections
            person  = unidrnd(N);
            while ( ~(isempty(find(PopWhoToInfect(index,:) == person))) || person == index )
                person = unidrnd(N);
            end

            PopWhoToInfect(index, connectM ) = unidrnd(N);


        end
    end
    
end



% establishing the network beforehand and keeping it
% 
% If someone has had this that or the other virus, then they can not get
% the newly mutated virus
% Virus bit string 0s to 1s and 1s to 0s only
% Change a random bit in the bit string, use unidrand(7)

numberOfInfected = zeros(tend,1);
PopCanInfect = zeros(N,1);
tout = 1:tend;


for ( t = 1:tend )
   
    t
    
    infectedCount = numel(find(Pop == 1));
    
        
    if ( infectedCount == 0) %this stops the main loop if no people are infected, as 
        for ( r = t : tend)  %the virus is entirely dead and there is no point in continuing the loop
            
            numberOfInfected(r,1) = 0;
        end
        break;

    end
    
    
    
    numberOfInfected(t,1) = infectedCount;
    
    
    
    
    for ( n = 1:N)
        if ( Pop(n) == 1)             % Trying to make the virus only move from an 
            PopCanInfect(n) = 1;      % infected to a susceptible once per time step 
                                      % (person A infects B and then B can't infect until t = t + 1 )
        end
    end
    
     for ( n = 1:N)  %Infecting everyone near an infected person, with a probability of infectProb
        if ( Pop(n) == 1 && PopCanInfect(n) == 1)  %making sure that only people who were infected at the beginning
            for i = 1 : 2 * K                      %of the time step infect people during this time step
                if (rand(1) < infectProb && Pop(PopWhoToInfect(n,i)) == 0 ) % This assumes that a person who is
                        %can a person who is currently infected with a
                        %virus, can NOT get infected with a NEW virus (even
                        %if they are not immune)
                    if ( isempty(VH(PopWhoToInfect(n,i)).viruses) == 1 ) %if the person has had no viruses previously
                        
                        if ( rand(1) < mu ) %probability of the virus mutating, the virus would
                            spot = unidrnd(virusSize); %change one bit inside the person who is transmitting it
                            if ( VH(n).virus(1, spot ) == 1 )
                                VH(n).virus(1,spot) = 0;
                            end
                            
                            if ( VH(n).virus(1, spot ) == 0 )
                                VH(n).virus(1,spot) = 1;
                            end
                            
                        end
                        
                        
                        %this moves the infection from person n to person 
                        % PopWhoToInfect(n,i)
                        Pop(PopWhoToInfect(n,i)) = 1;
                        VH(PopWhoToInfect(n,i)).viruses = VH(n).virus(1,:);
                        VH(PopWhoToInfect(n,i)).time = 0;
                        VH(PopWhoToInfect(n,i)).times =  0 ;
                        VH(PopWhoToInfect(n,i)).virus =  VH(n).virus(1,:);

                    end
                    
                    %if the person has any viruses in their history, then
                    %the hamming distances need to be checked
                    if ( ~( isempty(VH(PopWhoToInfect(n,i)).viruses) == 1) )
                        
                        
                        
                        num = size(VH(PopWhoToInfect(n,i)).viruses);
                        
                        hammingDistances = num(2) * ones(1,num(1)); %finding the hamming distance

                        for ( x = 1: num(1) )
                            for ( y = 1:num(2))
                                
                                if ( VH(PopWhoToInfect(n,i)).viruses(x,y) == VH(n).virus(1,y) )
                                    hammingDistances(1,x) = hammingDistances(1,x) - 1;

                                end
                            end
                        end

                        hammingDistance = min(hammingDistances);

                        if ( hammingDistance > hthr ) %infecting the person if they are not immune, that is, if the
                            
                            if ( rand(1) < mu )
                            spot = unidrnd(virusSize);
                                if ( VH(n).virus(1, spot ) == 1 )
                                    VH(n).virus(1,spot) = 0;
                                end
                            
                                if ( VH(n).virus(1, spot ) == 0 )
                                    VH(n).virus(1,spot) = 1;
                                end
                            
                            end
                            
                            Pop(PopWhoToInfect(n,i)) = 1; %viruses in their history are different enough from the current virus
                            VH(PopWhoToInfect(n,i)).viruses = [VH(PopWhoToInfect(n,i)).viruses ; VH(n).virus(1,:)];
                            VH(PopWhoToInfect(n,i)).time = 0;
                            VH(PopWhoToInfect(n,i)).times = [VH(PopWhoToInfect(n,i)).times ; 0];
                            VH(PopWhoToInfect(n,i)).virus = VH(n).virus(1,:);
                                                    
                        end                         
                    end   
                end
            end
        end
    end      

    
    
    for ( n = 1 : N)  %Increasing the time the infected people have been infected
          
        if (Pop(n) == 1 & VH(n).time > tauI)
            Pop(n) = 0;
            VH(n).virus = [];
            VH(n).time = [];
        
        elseif ( Pop(n) == 1)
            
            
            VH(n).time = VH(n).time + 1;
            
           
        end
      
        
    end
    

    for ( n = 1 : N)  %Increasing the time the people have had each virus they have had
            
        if ( isempty( find( VH(n).times > tauR)) == 0 ) 
            endTimes = find( VH(n).times > tauR );
            
            VH(n).viruses(endTimes,:) = [];
            VH(n).times(endTimes,:) = [];
        end
            howBig = size(VH(n).times);
            
            if ( howBig(2) == 0 )
                howBig(1) = howBig(2);
            end
        
        
        for ( k = 1 : howBig(1) )
                VH(n).times(k) =  VH(n).times(k) + 1;
        end
        
        
    end

    
    
    PopCanInfect = zeros(N,1);

    
end

toc

time = toc;

end
