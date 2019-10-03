function refpoint = imagesc4D( Img, disprange, f_ginput, initstr)
%imagesc4D displays 4D grayscale or RGB images in slice by slice fashion
%with mouse based slice browsing and window and level adjustment control.
%
% Usage:
% imagesc4D ( Image )
% imagesc4D ( Image , [] )
% imagesc4D ( Image , [LOW HIGH] )
% imagesc4D ( Image , 'ginput' )
% imagesc4D ( Image , 'ginput' , initstr)
% imagesc4D ( Image , [LOW HIGH], 'ginput' )
% imagesc4D ( Image , [LOW HIGH], 'ginput' , initstr)
%   
%   Image:      4D image MxNxKxLxC (L 3D volumes of K slices of MxN images) 
%               C is either 1 (for grayscale images) or 3 (for RGB images)  
%   [LOW HIGH]: display range that controls the display intensity range of
%               a grayscale image (default: the widest available range)
%   'ginut':    holds the function execution until 'S' is pressed on the
%               keyboard. When pressed, imagesc4D returns the 4 coordinates
%               of the voxel under the mouse pointer.
%   initstr:    to use with 'ginut'. Should contain all figure 
%               parametrization commands, such as colormap, axis, title, to
%               be execuded before the function imagesc4D waits for user
%               input.
%
% Use the scroll bar or mouse scroll wheel to switch between slices. Use 
% the second scroll bar or SHIFT + mouse scroll to switch between slabs 
% in the fourth dimension. To adjust window and level values keep the mouse 
% right button pressed and drag the mouse up and down 
% (for level adjustment) or right and left (for window adjustment). Window 
% and level adjustment control works only for grayscale images.
% 
% "Auto W/L" button adjust the window and level automatically for grayscale
% images
%
% While "Fine Tune" check box is checked the window/level adjustment gets
% 16 times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse based window and level adjustment is set
% based on the user defined display intensity range; the wider the range
% the more sensitivity to mouse drag.
% 
% 
%   Example
%   --------
%       % Display a 4D image
%       [X,Y,Z,T] = ndgrid(0:299,0:299,0:10,0:10);
%       im_phase4d = 2*pi*(X*0.2 + Y*0.1 + Z*1 + T*2);
%       figure();imagesc4D(im_phase4d);colormap(gray);
%
%       % Display the image, adjust the display range
%       figure();imagesc4D(im_phase4d,[-10 10]);colormap(gray);
%

%
% - Maysam Shahedi (mshahedi@gmail.com) (IMSHOW3D)
% - Anton Abyzov (ivoreus@gmail.com) (imagesc4D)
% - Released: 1.0.0   Date: 2018/06/01
% 

ExitVal = 0;

sno = size(Img,3);  % number of slices
eno = size(Img,4);  % number of echoes
S = round(sno/2);
E = 1;
im_sz = size(Img);


global InitialCoord;

MousePosition = [round(im_sz(2)/2) round(im_sz(1)/2)];

MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end    

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
else
    if (ischar(disprange) == 1)
        if exist('f_ginput')
            initstr = f_ginput;
        end
        f_ginput = disprange;
        disprange = [];
    end
    if numel(disprange) == 0
        [Rmin Rmax] = WL2R(Win, LevV);
    else
        LevV = (double(disprange(2)) + double(disprange(1))) / 2;
        Win = double(disprange(2)) - double(disprange(1));
        WLAdjCoe = (Win + 1)/1024;
        [Rmin Rmax] = WL2R(Win, LevV);
    end
end

if exist('f_ginput')
    if strcmp(f_ginput,'ginput')
        f_ginput = 1;
    else
        f_ginput = 0;
    end
else
    f_ginput = 0;
end
    
scrsz = get(0,'ScreenSize');
fig = gcf;
fig.Position = [scrsz(3)/3 0.15*scrsz(4) scrsz(3)/3 0.7*scrsz(4)];

h=axes('position',[0.1,0.25,0.85,0.7]);
imagesc(squeeze(Img(:,:,S,E,:)), [Rmin Rmax]); axis image;

FigPos = get(gcf,'Position');
if sno == 1
    epos2 = 45;
else
    epos2 = 85;
end
E_Pos = [50 epos2 uint16(FigPos(3)-100)+1 20];
Etxt_Pos = [50 epos2+20 uint16(FigPos(3)-100)+1 15];
S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [50 20 60 20];
Wval_Pos = [110 20 45 20];
Ltxt_Pos = [160 20 45 20];
Lval_Pos = [205 20 45 20];
BtnStPnt = uint16(FigPos(3)-250)+1;
if BtnStPnt < 300
    BtnStPnt = 300;
end
Sel_Pos = [BtnStPnt-40 20 60 20];
Btn_Pos = [BtnStPnt+30 20 80 20];
ChBx_Pos = [BtnStPnt+120 20 80 20];

if ((sno == 1) & (eno == 1))
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end
if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end
if eno > 1
    ehand = uicontrol('Style', 'slider','Min',1,'Max',eno,'Value',E,'SliderStep',[1/(eno-1) 10/(eno-1)],'Position', E_Pos,'Callback', {@SliceSlider2, Img});
    etxthand = uicontrol('Style', 'text','Position', Etxt_Pos,'String',sprintf('4D Slab# %d / %d',E, eno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
if f_ginput == 1
    Selhand = uicontrol('Style', 'pushbutton','Position', Sel_Pos,'String','Select..', 'FontSize', BtnSz, 'Callback' , @SelectPos);
end
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)
set (gcf,'WindowKeyPressFcn', @keyPress);
set (gcf,'WindowKeyReleaseFcn', @keyRelease);
set(gcf, 'WindowButtonMotionFcn', @mouseMoveLevel);

if exist('initstr')
    eval(initstr);
end

    function keyPress(object, eventdata)
        if strcmp(eventdata.Key,'shift')
            set (gcf, 'WindowScrollWheelFcn', @mouseScroll2);
        end
        if strcmp(eventdata.Key,'s')
            rowref = floor(MousePosition(2));
            colref = floor(MousePosition(1));
            refpoint = [rowref,colref,S,E];
            ExitVal = 1;
        end
    end

    function keyRelease(object, eventdata)
        if strcmp(eventdata.Key,'shift')
            set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
        end
    end

% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        E_Pos = [50 epos2 uint16(FigPos(3)-100)+1 20];
        Etxt_Pos = [50 epos2+20 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-250)+1;
        if BtnStPnt < 300
            BtnStPnt = 300;
        end
        Sel_Pos = [BtnStPnt-40 20 60 20];
        Btn_Pos = [BtnStPnt+30 20 80 20];
        ChBx_Pos = [BtnStPnt+120 20 80 20];
        if ((sno == 1) & (eno == 1))
            set(stxthand,'Position', Stxt_Pos);
        end
        if sno > 1
            set(shand,'Position', S_Pos);
            set(stxthand,'Position', Stxt_Pos);            
        end
        if eno > 1
            set(ehand,'Position', E_Pos);
            set(etxthand,'Position', Etxt_Pos);
        end
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        if f_ginput == 1
            set(Selhand,'Position', Sel_Pos);
        end
        set(ChBxhand,'Position', ChBx_Pos);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,E,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        end
        if ((sno == 1) & (eno == 1))
            set(stxthand, 'String', '2D image');
        end
    end

    function SliceSlider2 (hObj,event, Img)
        E = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,E,:)))
        caxis([Rmin Rmax])
        if eno > 1
            set(etxthand, 'String', sprintf('4D Slab# %d / %d',E, eno));
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        end
        if ((sno == 1) & (eno == 1))
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,E,:)))
    end

    function mouseScroll2 (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        E = E - UPDN;
        if (E < 1)
            E = 1;
        elseif (E > eno)
            E = eno;
        end
        if eno > 1
            set(ehand,'Value',E);
            set(etxthand, 'String', sprintf('4D Slab# %d / %d',E, eno));
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,E,:)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', @mouseMoveLevel);
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        try
            cursorPoint = get(gca, 'CurrentPoint');
            MousePosition = [cursorPoint(1,1) cursorPoint(1,2)];
        end
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

    function mouseMoveLevel(object, eventdata, handles)
        try
            cursorPoint = get(gca, 'CurrentPoint');
            MousePosition = [cursorPoint(1,1) cursorPoint(1,2)];
        catch
            %disp('Mouse moved outside bounds');
        end
    end


% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        try
            cursorPoint = get(gca, 'CurrentPoint');
            MousePosition = [cursorPoint(1,1) cursorPoint(1,2)];
        end
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

    function SelectPos(object,eventdata)
        [colref,rowref] = ginput(1); 
        colref = round(colref); rowref = round(rowref);
        refpoint = [rowref,colref,S,E];
        ExitVal = 1;
    end

    if (f_ginput == 1)
        while (ExitVal == 0)
            drawnow;
        end
        disp('Exit');
    end

end
% -=< Maysam Shahedi (mshahedi@gmail.com), September 22, 2016>=-
% -=< Anton Abyzov (ivoreus@gmail.com), June 01, 2018>=-