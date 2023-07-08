function saveFigs(targetdir)

a=dir(strcat(pwd,'/*.fig'));
for i=1:length(a),
    c=strsplit(a(i).name,'.');
    fig=openfig(a(i).name);
    f=sprintf('%s/%s.eps',targetdir,c{1});
    print(fig,f,'-depsc2');
    close(fig);
end