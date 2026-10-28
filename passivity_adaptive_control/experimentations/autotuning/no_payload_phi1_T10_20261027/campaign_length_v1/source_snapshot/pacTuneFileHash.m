function h=pacTuneFileHash(file)
%PACTUNEFILEHASH Streaming SHA-256; original file is opened read-only.
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(file,'rb');assert(fid>=0,'pacTune:File','Cannot open %s.',file);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
while ~feof(fid),md.update(fread(fid,1048576,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
