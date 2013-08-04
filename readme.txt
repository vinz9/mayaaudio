------------------------------------------------
install :

copy the version of the plugin corresponding to your maya version (right now it is compiled for Maya 2009, x86 and x64 versions) in a folder called plug-ins in your maya documents folder.
copy libfftw3f-3.dll in the bin directory of your maya installation (where maya.exe is located)


------------------------------------------------
usage :

check the mel procedures to get you started.
only 8 and 16 bits uncompressed wav files, mono or stereo, are supported, since it is what is supported by Maya.




------------------------------------------------
parameters :

Bands Number : number of bands the audio is split into
Amp Scale : scales the amplitude
Amp response : speed at which whatever you plug in the outputs reacts (values between 0.9 and 0.99 are usually good). If everything goes weird, putting it back to 0 temporarily should fix problems.
Fft scaling : various ways to scale the spectrum to better reflect human perception
Samples Number : precision of the fourier Transform in samples


-------------------------------------------------

You can change the file path and the offset on the maya audio node to which the audioNode node (sorry, bad naming convention) is connected.
You should bake the animation and unload the plugin when you're done.
Not much error checking is done, so save your file before using the plugin, it can crash with bad audio files.



04/08/2013
12/09/2012 : version 0.91 : compiled for Maya 2011/2012/2013
19/04/2009 : version 0.9 : initial (and probably last) public release, after more than one year sitting on my hard drive

Please credit me and drop me a line if you do something cool with it!

vincent.houze@foliativ.net
http://www.foliativ.net



