using Distributions
using FreqTables
using Gadfly
using Cairo

#This Julia script can be used to generate random (binomially)
#distributed codon usage bias data, which can then be compared to
#the actual sequences.
#It contains a function (getSample) to produce the random data.
#The cub can be adjusted in the function body
#So can the length and the noise around the cub.
#The script also contains a few lines to evaluate the slope

#############################################
#Calculates artificial (random) Entropy data for a 
#probability "bas"
#Returns the data that can then be used to calulcated the 
#temperature
#If basenorm is true, then entropy is calculated
#relative to equal opportunity model
#############################################
function getSample(bas=0.5,basenorm=true)
#tmp is the success probability
#rt is the length we consider
local tmp5,tmp3,tmp4,tbl,outp,v,d,based
 tmp5,tmp3,tmp4,tbl,outp=-1,-1,-1,-1,-1

v=0;
for i in 1:1

#first decide on the size N-> rt and the probablity p->tmp 
rt = rand([5,6,7,8,9,10])
rt=5
#tmp = 0.5 + 0.1*randn()
tmp =bas + 0.*randn()

if (tmp > 1) tmp =1 end
if(tmp <0) tmp = 0 end



#Now define your distributions d is the biased one, based is the one with 0.5
d=Binomial(rt,tmp)

#if basenorm is true, then use as base probability the 0.5 prob otherwise the same as the input prob
if(basenorm)
based= Binomial(rt,0.5)
else
based= Binomial(rt,tmp)
end

#DELETEME DELTE NEXT LINE
#based= Binomial(rt,0.2)


#inner loop
#generate samples with the chosen probability and size
	
for j in 1:1

#generate the actual sample
tmp10=rand(d,10000)

#calculate the probability
a=round.(1000000000*pdf.(based,tmp10))/1000000000


if ( ((rt % 2) == 0) & (tmp10 ==  Int64(round(rt/2))  ))
#println("Twice")
v=cat(1,v,a)
else

if(basenorm)
v=cat(1,v,2*a)
else # don't mulyoply by 2
v=cat(1,v,a)
end#if

end#if

end#innerfor

#now produce the freqtable
#This freqtable has the frequency as col2
# and the probability as col1
deleteat!(v,1)
a= v
tbl=freqtable(a)
tmp4=convert(Array, tbl)
en=sum(tmp4)
rel=names(tbl)[1]
tmp4=tmp4/sum(tmp4)# sampled frequencies (empirical distribution)
#tmp6= log.(tmp4)
tmp6= tmp4
#tmp3=-log.(names(tbl)[1])
tmp3=names(tbl)[1]#theoretical distributions
#println(tbl)
println(tmp4)
println( exp(-en*sum(tmp4.*log.(tmp4./rel)) ))
println( sum(tmp4.*log.(tmp4./rel)) )
tmp5=[tmp3 tmp6]
if(outp != -1)
outp = [outp; tmp5]
else
outp = tmp5
end

v=0
#Now start the whole thing again, by choosing a new number and length
end#outerfor


return(outp)
end#function
#############################################
#############################################

#if you pass false to this function, then the probability is evaluated 
#according to the same bias as the sample sequence.
#Otheriwise to the binomial with p=1/2
a= getSample(0.2,false)

#writedlm("/tmp/dta.txt",a)

writecsv("/tmp/dta.txt",a)
tmp=linreg(a[:,1],a[:,2])
print(tmp)
#myplot=plot(x=a[:,1],y=a[:,2],Geom.point,Geom.smooth(method=:lm))
#draw(PS("myplot.eps", 8inch, 6inch), myplot)

