function obj = recordOnFilm(obj)
% Record the psf onto the film from the rays
%
% Uses the current ppsfRays and the film objects.  These are
% both part of the camera object.
%
% Example:
%    obj.rays.recordOnFilm(film)

obj.rays.recordOnFilm(obj.film);

end
