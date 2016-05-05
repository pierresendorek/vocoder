using WAV

x,Fs= wavread("in.wav")

fqA4=440.0

type Phase
  phi::Float64 # dérivée de fréquence instantanée en Hz
end

function freqToNote(f::Float64)
  return log(f/fqA4)*12/log(2)
end

function noteToFreq(x::Float64)
  return fqA4*exp(x*log(2)/12)
end

function incrementPhaseByFreq!(phase::Phase,freq::Float64)
  phase.phi = rem( phase.phi + freq/Fs , 1.0 )
end

function saw(phase::Phase)
  return 2*(phase.phi-0.5)
end

y = zeros(Float64,length(x))

noteArray=[-24.0,-24.0+4.0,-27.0+7.0,0.0,4.0,7.0]

phaseArray=Array(Phase,length(noteArray))
for iPhase in 1:length(phaseArray)
  phaseArray[iPhase]=Phase(0)
end

for it in 1:length(y)
  for iNote in 1:length(noteArray)
    y[it]+= saw(phaseArray[iNote])
    incrementPhaseByFreq!(phaseArray[iNote],noteToFreq(noteArray[iNote]))
  end
end

y=y/maximum(abs(y))
wavwrite(y,"carrier.wav",Fs=Fs)


z=zeros(Float64,length(x))

winLen = 1024
apoWin = sinpi(linspace(0,1,winLen)).^2
gap=512


iWin=0
sigma2=0.001

while iWin*gap+1 + winLen <= length(x)
  X=fft(x[iWin*gap+1:iWin*gap+winLen].*apoWin)
  Y=fft(y[iWin*gap+1:iWin*gap+winLen].*apoWin)
  tY=Y.*abs(X).*abs(Y)./(abs(Y).^2 + sigma2)
  z[iWin*gap+1:iWin*gap+winLen]+=real(ifft(tY))
  iWin+=1
end

z=z/maximum(abs(z))

wavwrite(z,"out.wav",Fs=Fs)







