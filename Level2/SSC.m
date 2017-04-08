function [ frameType ] = SSC( frameT, nextFrameT, prevFrameType )
%FrameType: OLS, LSS, ESH, LPS
%frameT: Current frame - 2048x2 matrix
%nextFrameT: next frame - 2048x2 matrix
%prevFrameType: The type of the previous frame

%Apply filter on next frame
b = 0.7548 * [1, -1];
a = [1, -0.5095];

y = filter(b, a, nextFrameT);

%Create non-overlapping frames of length 128 samples
%frames is 128x16 matrix - first 8 columns are the channel 1 subframes
%and columns 9-16 the channel 2 subframes

frames = [buffer(y(577:end-448, 1), 128) ...
          buffer(y(577:end-448, 2), 128)];

%energy is a 8x2 matrix - 1 column per channel - 8 subframes
%same for ds but here we have 7 subframes
energy = reshape(sum(frames.^2), 8, 2);

temp = cumsum(energy);
ds =  [1:7; 1:7]' .* energy(2:end, :) ./ temp(1:end-1, :);

eight_logic = ( energy(2:end, :) > 1e-3 ) & ( ds > 10 );

%nextFrameType = TRUE if ESH else FALSE
% nextFrameType1 = any(eight_logic(:, 1));
% nextFrameType2 = any(eight_logic(:, 2));
% 
% types = ['OLS' 'LSS' 'ESH' 'LSS' 'LSS' 'LSS' 'ESH' 'ESH' 'ESH' 'ESH' 'ESH' 'ESH' 'LPS' 'ESH' 'ESH' 'LPS'];

%Decide if the next frame is ESh or not
next_is_eight = any(any(eight_logic));

if (strcmp(prevFrameType, 'OLS') && next_is_eight)
    frameType = 'LSS';
elseif (strcmp(prevFrameType, 'OLS') && (~next_is_eight))
    frameType = 'OLS';
elseif (strcmp(prevFrameType, 'ESH') && next_is_eight)
    frameType = 'ESH';
elseif (strcmp(prevFrameType, 'ESH') && (~next_is_eight))
    frameType = 'LPS';
elseif (strcmp(prevFrameType, 'LSS'))
    frameType = 'ESH';
elseif (strcmp(prevFrameType, 'LPS'))
    frameType = 'OLS';
end

% if (strcmp(prevFrameType, 'OLS') && nextFrameType1)
%     frameType1 = 'LSS';
% elseif (strcmp(prevFrameType, 'OLS') && (~nextFrameType1))
%     frameType1 = 'OLS';
% elseif (strcmp(prevFrameType, 'ESH') && nextFrameType1)
%     frameType1 = 'ESH';
% elseif (strcmp(prevFrameType, 'ESH') && (~nextFrameType1))
%     frameType1 = 'LPS';
% elseif (strcmp(prevFrameType, 'LSS'))
%     frameType1 = 'ESH';
% elseif (strcmp(prevFrameType, 'LPS'))
%     frameType1 = 'OLS';
% end
% 
% if (strcmp(prevFrameType, 'OLS') && nextFrameType2)
%     frameType2 = 'LSS';
% elseif (strcmp(prevFrameType, 'OLS') && (~nextFrameType2))
%     frameType2 = 'OLS';
% elseif (strcmp(prevFrameType, 'ESH') && nextFrameType2)
%     frameType2 = 'ESH';
% elseif (strcmp(prevFrameType, 'ESH') && (~nextFrameType2))
%     frameType2 = 'LPS';
% elseif (strcmp(prevFrameType, 'LSS'))
%     frameType2 = 'ESH';
% elseif (strcmp(prevFrameType, 'LPS'))
%     frameType2 = 'OLS';
% end
% 
% keySet = {'OLS', 'LSS', 'ESH', 'LPS'};
% map1 = containers.Map(keySet, 0:3);
% map2 = containers.Map(keySet, [1 4 7 10]);
% 
% temp = map1(frameType1) * 12 + map2(frameType2);
% 
% frameType = types(temp:temp+2);
% 
% if (strcmp(frameType1, 'ESH') || strcmp(frameType2, 'ESH'))
%     frameType = 'ESH';
% elseif (strcmp(frameType1, 'LSS') && strcmp(frameType2, 'LPS')) || (strcmp(frameType1, 'LPS') && strcmp(frameType2, 'LSS'))
%     frameType = 'ESH';
% elseif (strcmp(frameType1, 'LSS') || strcmp(frameType2, 'LSS'))
%     frameType = 'LSS';
% elseif (strcmp(frameType1, 'LPS') || strcmp(frameType2, 'LPS'))
%     frameType = 'LPS';
% else
%     frameType = 'OLS';
% end

end

