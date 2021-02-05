/* global ol, view, map */

var view = new ol.View({
    center: ol.proj.transform([20.9, 44.0], 'EPSG:4326', 'EPSG:3857'),
    zoom: 7});
(function () {
    var container = document.getElementById('popup');
    var content = document.getElementById('popup-content');
    var closer = document.getElementById('popup-closer');
    var overlay = new ol.Overlay(/** @type {olx.OverlayOptions} */ ({
        element: container,
        autoPan: true,
        autoPanAnimation: {
            duration: 250
        }
    }));
    closer.onclick = function () {
        overlay.setPosition(undefined);
        closer.blur();
        return false;
    };
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
            new ol.layer.Group({
                title: 'Podloga',
                layers: [
                    new ol.layer.Group({
                        title: 'Water Color sa oznakama',
                        type: 'base',
                        combine: true,
                        visible: false,
                        layers: [
                            new ol.layer.Tile({
                                source: new ol.source.Stamen({
                                    layer: 'watercolor'
                                })
                            }),
                            new ol.layer.Tile({
                                source: new ol.source.Stamen({
                                    layer: 'terrain-labels'
                                })
                            })
                        ]
                    }),
                    new ol.layer.Tile({
                        title: 'Natural Earth',
                        type: 'base',
                        combine: true,
                        visible: false,
                        source: new ol.source.TileJSON({
                            url: 'https://api.tiles.mapbox.com/v3/mapbox.natural-earth-hypso-bathy.json?secure'
                        })
                    }),
                    new ol.layer.Tile({
                        title: 'Open Street Map',
                        type: 'base',
                        visible: true,
                        source: new ol.source.OSM()
                    })
                ]
            }),
            new ol.layer.Group({
                title: 'Lejeri',
                layers: [
                    new ol.layer.Tile({
                        title: 'Erozija (Geoserver)',
                        visible: false,
                        source: new ol.source.TileWMS({
                            url: 'http://localhost:8080/geoserver/Mirko/wms',
                            params: {'LAYERS': 'Mirko:erozija'},
                            serverType: 'geoserver'
                        })
                    }),
                    new ol.layer.Tile({
                        title: 'Stanice (Geoserver)',
                        visible: false,
                        source: new ol.source.TileWMS({
                            url: 'http://localhost:8080/geoserver/Mirko/wms',
                            params: {'STYLES': 'stanice', 'LAYERS': 'Mirko:stanice'},
                            serverType: 'geoserver'
                        })
                    })
                ]
            })
        ],
        overlays: [overlay],
        view: view
    });
    var layerSwitcher = new ol.control.LayerSwitcher({
        tipLabel: 'Légende' // Optional label for button
    });
    map.addControl(layerSwitcher);
    var highlightStyleCache = {};
    var featureOverlay = new ol.layer.Vector({
        source: new ol.source.Vector(),
        map: map,
        style: function (feature, resolution) {
            var text = resolution < 5000 ? feature.get('name') : '';
            if (!highlightStyleCache[text]) {
                highlightStyleCache[text] = new ol.style.Style({
                    stroke: new ol.style.Stroke({
                        color: '#f00',
                        width: 1
                    }),
                    fill: new ol.style.Fill({
                        color: 'rgba(255,0,0,0.1)'
                    }),
                    text: new ol.style.Text({
                        font: '12px Calibri,sans-serif',
                        text: text,
                        fill: new ol.style.Fill({
                            color: '#000'
                        }),
                        stroke: new ol.style.Stroke({
                            color: '#f00',
                            width: 3
                        })
                    })
                });
            }
            return highlightStyleCache[text];
        }
    });
    var highlight;
    var displayFeatureInfo = function (pixel) {
        var feature = map.forEachFeatureAtPixel(pixel, function (feature) {
            return feature;
        });
        var info = document.getElementById('info');
        if (feature) {
            info.innerHTML = feature.getId() + ': ' + feature.get('name');
        } else {
            info.innerHTML = '&nbsp;';
        }
        if (feature !== highlight) {
            if (highlight) {
                featureOverlay.getSource().removeFeature(highlight);
            }
            if (feature) {
                featureOverlay.getSource().addFeature(feature);
            }
            highlight = feature;
        }
    };
    map.on('pointermove', function (evt) {
        if (evt.dragging) {
            return;
        }
        var pixel = map.getEventPixel(evt.originalEvent);
        displayFeatureInfo(pixel);
    });
    map.on('click', function (evt) {
        displayFeatureInfo(evt.pixel);
    });
    var style = new ol.style.Style({
        fill: new ol.style.Fill({
            color: 'rgba(255, 255, 255, 0.6)'
        }),
        stroke: new ol.style.Stroke({
            color: '#319FD3',
            width: 1
        }),
        text: new ol.style.Text({
            font: '12px Calibri,sans-serif',
            fill: new ol.style.Fill({
                color: '#000'
            }),
            stroke: new ol.style.Stroke({
                color: '#fff',
                width: 3
            })
        })
    });
    map.on('singleclick', function (evt) {
        var coordinate = evt.coordinate;
        var hdms = ol.coordinate.toStringHDMS(ol.proj.transform(
                coordinate, 'EPSG:3857', 'EPSG:4326'));
        content.innerHTML = '<p>Koordinate pozicije:</p><code>' + hdms +
                '</code>';
        overlay.setPosition(coordinate);
    });
    function bindInputs(layerid, layer) {
        var visibilityInput = $(layerid + ' input.visible');
        visibilityInput.on('change', function () {
            layer.setVisible(this.checked);
        });
        visibilityInput.prop('checked', layer.getVisible());
        var opacityInput = $(layerid + ' input.opacity');
        opacityInput.on('input change', function () {
            layer.setOpacity(parseFloat(this.value));
        });
        opacityInput.val(String(layer.getOpacity()));
    }

    map.getLayers().forEach(function (layer, i) {
        bindInputs('#layer' + i, layer);
        if (layer instanceof ol.layer.Group) {
            layer.getLayers().forEach(function (sublayer, j) {
                bindInputs('#layer' + i + j, sublayer);
            });
        }
    });
    $('#layertree li > span').click(function () {
        $(this).siblings('fieldset').toggle();
    }).siblings('fieldset').hide();



    //CRTANJE

    // init a vector layer whose source can be used with our new control
    var vectorSource = new ol.source.Vector();

    var vectorLayer = new ol.layer.Vector({
        source: vectorSource,
        style: new ol.style.Style({
            fill: new ol.style.Fill({
                color: 'rgba(255, 255, 0, 0.3)'
            }),
            stroke: new ol.style.Stroke({
                color: '#ffff00',
                width: 3
            }),
            image: new ol.style.Circle({
                radius: 5,
                fill: new ol.style.Fill({
                    color: '#ffff00'
                })
            })
        })
    });

    map.addLayer(vectorLayer); // here we use our new control

    var typeSelect = document.getElementById('type');

    var draw; // global so we can remove it later

    function addInteraction() {
        var value = typeSelect.value;
        if (value !== 'None') {
            draw = new ol.interaction.Draw({
                source: vectorSource,
                type: /** @type {ol.geom.GeometryType} */ (typeSelect.value)
            });
            map.addInteraction(draw);
        }
    }

    /* Handle change event. */
    typeSelect.onchange = function () {
        map.removeInteraction(draw);
        addInteraction();
    };

    addInteraction();
})();



