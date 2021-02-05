/* global ol, google, view */


var map = new ol.Map({
    target: 'map',
    renderer: 'canvas',
    controls: ol.control.defaults({
        attributionOptions: ({
            collapsible: false
        })
    }).extend([
        new ol.control.ZoomSlider(),
        new ol.control.ZoomToExtent({}),
        new ol.control.Rotate(),
        new ol.control.OverviewMap(),
        new ol.control.ScaleLine(),
        new ol.control.MousePosition({
            coordinateFormat: ol.coordinate.createStringXY(4),
            projection: 'EPSG:4326'
        })
    ]),
    // interactions and controls are seperate entities in ol3
    // we extend the default navigation with a hover select interaction
    interactions: ol.interaction.defaults().extend([
        new ol.interaction.Select({
            condition: ol.events.condition.mouseMove
        })
    ]),
    layers: [
        new ol.layer.Tile({
            source: new ol.source.OSM()
        })
    ],
    view: view
});



// Pretraga lokacije, zoom level

function lonlat() {
    var lat = document.getElementById("lat").value;
    var lon = document.getElementById("lon").value;
    var center = [parseInt(lon), parseInt(lat)];
    map.getView().setCenter(ol.proj.transform(center, 'EPSG:4326', 'EPSG:3857'));
}

function GetLocation() {
    var geocoder = new google.maps.Geocoder();
    var address = document.getElementById("adresa").value;
    geocoder.geocode({'address': address}, function (results, status) {

        var latitude = results[0].geometry.location.lat();
        var longitude = results[0].geometry.location.lng();

        var coord = [longitude, latitude];
        var centre = ol.proj.transform(coord, 'EPSG:4326', 'EPSG:3857');
        var map2 = view.setCenter(centre);
        var zoom2 = view.setZoom(17);

    });

}

function setLevel() {
    var level = document.getElementById("level").value;
    map.getView().setZoom(level).value;
}
