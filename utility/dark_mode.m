function dark_mode(use_dark_mode)
  % DARK_MODE  Switch command and figures windows to dark background with light
  % text.
  %
  % dark_mode
  % dark_mode(flag)
  %
  % Optional input:
  %   flag  whether to switch to dark mode {false} or back to default (true)
  %
  % Examples:
  %   % Place this in startup.m to add a timer that checks every 5mins whether
  %   % its ≥6am or ≥6pm and toggles dark mode accordingly.
  %   dmt = timer('TimerFcn',@(x,y) dark_mode(mod(clock*[0;0;0;1;0;0]+6,24)<12),'Period',60*5,'ExecutionMode','fixedDelay');
  %   start(dmt);
  %
  %   % Detect if Mac OS X is dark and match
  %   dark_mode(system('defaults read -g AppleInterfaceStyle 2&>/dev/null')==0);
  %

  % Adapted from:
  % https://github.com/benhager/solarized-matlab/blob/master/setupSolarized.m

  if ~exist('use_dark_mode','var') || use_dark_mode
    sys = false;
    txc = [1 1 1];
    bgc = [30 30 30]/255;
    fbg = [1 1 1].*40/255;
    alc = [1 1 1];
    ac = [30 30 30]/255;
    str = [254.1188  106.3358  255.0000]/255;
    ustr = [203  75  22] / 255;
    hyp = [ 38 139 210] / 255;
    cmt = [ 54 139  54] / 255;
    scmd = [198 160  20] / 255;
    errs = [230 40 40]/255;
  else
    sys = true;
    txc = [0 0 0];
    bgc = [1 1 1];
    fbg = get(0,'factoryFigureColor');
    alc = 'k';
    ac = 'w';
    str = [160  32 240]/255;
    ustr = [178   0   0]/255;
    hyp = [0 0 1];
    cmt = [ 34 139  34] / 255;
    scmd = [178 140   0] / 255;
    errs = [230 0 0]/255;
  end
  kwd = hyp;
  warn = ustr;

  % Main colors
  com.mathworks.services.Prefs.setBooleanPref('ColorsUseSystem',sys);
  com.mathworks.services.Prefs.setColorPref('ColorsText',java.awt.Color(txc(1), txc(2), txc(3))); 
  com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsText');
  com.mathworks.services.Prefs.setColorPref('ColorsBackground',java.awt.Color(bgc(1), bgc(2), bgc(3))); 
  com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsBackground');
  com.mathworks.services.Prefs.setColorPref('Colors_M_Strings',java.awt.Color(str(1), str(2), str(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Strings');
  com.mathworks.services.Prefs.setColorPref('Colors_M_UnterminatedStrings',java.awt.Color(ustr(1), ustr(2), ustr(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_UnterminatedStrings');
  com.mathworks.services.Prefs.setColorPref('Colors_M_Keywords',java.awt.Color(kwd(1), kwd(2), kwd(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Keywords');
  com.mathworks.services.Prefs.setColorPref('Colors_M_Comments',java.awt.Color(cmt(1), cmt(2), cmt(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Comments');
  com.mathworks.services.Prefs.setColorPref('Colors_M_SystemCommands',java.awt.Color(scmd(1), scmd(2), scmd(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_SystemCommands');
  com.mathworks.services.Prefs.setColorPref('Colors_M_Errors',java.awt.Color(errs(1), errs(2), errs(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Errors');
  com.mathworks.services.Prefs.setColorPref('Colors_HTML_HTMLLinks',java.awt.Color(hyp(1), hyp(2), hyp(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_HTML_HTMLLinks');
  com.mathworks.services.Prefs.setColorPref('Colors_M_Warnings',java.awt.Color(warn(1), warn(2), warn(3)));
  com.mathworks.services.ColorPrefs.notifyColorListeners('Colors_M_Warnings');
  %com.mathworks.services.Prefs.setBooleanPref('ColorsUseMLintAutoFixBackground',afhb);
  %com.mathworks.services.Prefs.setColorPref('ColorsMLintAutoFixBackground',java.awt.Color(afh(1), afh(2), afh(3)));
  %com.mathworks.services.ColorPrefs.notifyColorListeners('ColorsMLintAutoFixBackground');
  %com.mathworks.services.Prefs.setColorPref('Editor.NonlocalVariableHighlighting.TextColor',java.awt.Color(vwss(1), vwss(2), vwss(3)));
  %com.mathworks.services.ColorPrefs.notifyColorListeners('Editor.NonlocalVariableHighlighting.TextColor');
  %com.mathworks.services.Prefs.setBooleanPref('EditorCodepadHighVisible',hsb);
  %com.mathworks.services.Prefs.setColorPref('Editorhighlight-lines', java.awt.Color(hsbc(1), hsbc(2), hsbc(3)));
  %com.mathworks.services.Prefs.setBooleanPref('Editorhighlight-caret-row-boolean',hclb);
  %com.mathworks.services.Prefs.setColorPref('Editorhighlight-caret-row-boolean-color',java.awt.Color(hcl(1), hcl(2), hcl(3)));
  %com.mathworks.services.ColorPrefs.notifyColorListeners('Editorhighlight-caret-row-boolean-color');
  %com.mathworks.services.Prefs.setBooleanPref('EditorShowLineNumbers',slnb);
  %com.mathworks.services.Prefs.setBooleanPref('EditorRightTextLineVisible',rtlb);

  % Figure windows
  set(0,'defaultfigurecolor',fbg)
  arrayfun(@(f) set(f,'Color',fbg),get(groot, 'Children'));
  % axes labels
  set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor','DefaultAxesColor','DefaultTextColor','DefaultLineColor'},{alc,alc,alc,ac,alc,alc})
  arrayfun(@(a) set(a,{'Xcolor','YColor','ZColor','Color'},{alc,alc,alc,ac}),findall(0,'type','axes'));
  arrayfun(@(t) set(t,'Color',alc),findall(0,'type','text'));



end
