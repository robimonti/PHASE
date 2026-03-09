function [municipality, country] = get_place_from_coordinates(lon, lat)

% get_place_from_coordinates Get municipality name and country for a given latitude and longitude

% API Endpoint URL
apiURL = sprintf('https://nominatim.openstreetmap.org/reverse?format=jsonv2&lat=%.6f&lon=%.6f', lat, lon);

% Web request options with timeout
options = weboptions('Timeout', 15);

try
% Make the web request
response = webread(apiURL, options);

% Extract municipality name (city, town, or village)
if isfield(response, 'address')
  if isfield(response.address, 'city')
    municipality = response.address.city;
  elseif isfield(response.address, 'town')
    municipality = response.address.town;
  elseif isfield(response.address, 'village')
    municipality = response.address.village;
  else
    municipality = 'Unknown';
  end
  if isfield(response.address, 'country')
    country = response.address.country;
  end
else
  municipality = 'Unknown (no address data in response)';
end
catch ME
% Handle potential errors during web request
disp('Error: Unable to retrieve place information.');
disp(ME.message);
municipality = 'Unknown';
country = 'Unknown';
end
end