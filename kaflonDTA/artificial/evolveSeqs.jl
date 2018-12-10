using Distributions
using FreqTables
using Gadfly
using Cairo

#a model of how *sequences* as opposed to
# individual codons evolve
#We implement a random walk model of
#the biomial distribution
#We a force can be added by biasing the probability
#to step to the left. see tag: ##BIAS
# for the relevant line of code.
# When the bias is 0, then the model just reproduces the binomial distribution

srand(1001)
seqLen=6
energs=Array{Float64}(seqLen+1)
energs[1]=0
rang=collect(0:seqLen)
seq=rand(rang)
seq=3
seq2= seqLen - seq
#sequence is (seq,seq2)
numtrial=100000# number of points we sample
relaxTime=2000# time we give the chain to relax
pow=2.0

#outer for
# for ty = collect(1)
rec=[]
for j = 1: numtrial


# we add some random elemnt to relax time to avoid artefacts 
numexp=relaxTime +Int64(round(rand()*100))# we add some random elemnt to relax time to avoid artefacts 

################
for i = 1:numexp
ppr= Float64((seqLen - seq)^pow) # we assume that seq is k. and on the left of the chain, k=0
ppl =Float64(seq)   

pl=ppl/(ppl+ppr)
pr=ppr/(ppl+ppr)

#seq2 is strictly for nothing, I include it anyway
#println(seq, " ", seq2, " ", pl, " " , pr)



if(rand() < pl)
seq = seq -1
seq2 = seq2 + 1
else
seq = seq +1
seq2=seq2 -1
end#if

#println(entropy, "  ", ppr, "   ",ppl, "  seq: ",seq )



end#inner for


#record the last value from last time
push!(rec,seq)
end#outer for


#Calculate the energies:
energs[1]=0;
for i = 2: seqLen+1
l=i-1
energs[i] = -log((seqLen-l +1)^pow/(l)) + energs[i-1]
end
println(energs)
Z=sum(exp.(-energs))
println(Z)
probs=(exp.(-energs))/Z
println(probs)





 outp= [collect(1:numtrial) rec]

#tbl= freqtable(Array{Float64,1}(round.(10000000*(pdf.(dis,rec)))/10000000))
tbl1= freqtable(Array{Float64,1}(rec))
#println(tbl1)
tmp=convert(Array, tbl1)
tmp=tmp/sum(tmp)
eprobs=tmp
println(tmp)
sequ=names(tbl1)[1]#sequences
cub=sum(tmp.*(sequ/seqLen))
#dis=Binomial(seqLen,cub)
dis=Binomial(seqLen,0.14322)
tprobs=pdf.(dis,names(tbl1)[1])


tmp=convert(Array, tbl1)
tmp=tmp/sum(tmp)
#println("Codon usage bias ",cub)

#plot(x=-tprobs,y=-eprobs,Geom.point,Geom.smooth(method=:lm))

tmp=linreg(-eprobs,-tprobs)
println(cub," ",tmp[2])
#end#outer for ty
