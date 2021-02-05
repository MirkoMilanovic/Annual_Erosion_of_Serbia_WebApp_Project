function applyMargins() {
    var leftToggler = $(".mini-submenu-left");
    var rightToggler = $(".mini-submenu-right");
    if (leftToggler.is(":visible")) {
        $("#map .ol-zoom")
                .css("margin-left", 0)
                .removeClass("zoom-top-opened-sidebar")
                .addClass("zoom-top-collapsed");
    } else {
        $("#map .ol-zoom")
                .css("margin-left", $(".sidebar-left").width())
                .addClass("zoom-top-opened-sidebar")
                .removeClass("zoom-top-collapsed");
    }

    if (leftToggler.is(":visible")) {
        $("#map .ol-zoomslider")
                .css("margin-left", 0)
                .removeClass("zoomslider-top-opened-sidebar")
                .addClass("zoomslider-top-collapsed");
    } else {
        $("#map .ol-zoomslider")
                .css("margin-left", $(".sidebar-left").width())
                .addClass("zoomslider-top-opened-sidebar")
                .removeClass("zoomslider-top-collapsed");
    }

    if (leftToggler.is(":visible")) {
        $("#map .ol-zoom-extent")
                .css("margin-left", 0)
                .removeClass("zoom-extent-top-opened-sidebar")
                .addClass("zoom-extent-top-collapsed");
    } else {
        $("#map .ol-zoom-extent")
                .css("margin-left", $(".sidebar-left").width())
                .addClass("zoom-extent-top-opened-sidebar")
                .removeClass("zoom-extent-top-collapsed");
    }

    if (leftToggler.is(":visible")) {
        $("#map .ol-overviewmap")
                .css("margin-left", 0)
                .removeClass("zoom-top-opened-sidebar")
                .addClass("zoom-top-collapsed");
    } else {
        $("#map .ol-overviewmap")
                .css("margin-left", $(".sidebar-left").width())
                .addClass("zoom-top-opened-sidebar")
                .removeClass("zoom-top-collapsed");
    }
    if (leftToggler.is(":visible")) {
        $("#map .ol-overviewmap")
                .css("margin-left", 0)
                .addClass("overviewmap-bottom");
    } else {
        $("#map .ol-overviewmap")
                .css("margin-left", $(".sidebar-left").width())
                .addClass("overviewmap-bottom");
    }
    if (leftToggler.is(":visible")) {
        $("#map .ol-scale-line")
                .css("margin-left", 0);
    } else {
        $("#map .ol-scale-line")
                .css("margin-left", $(".sidebar-left").width());
    }

    if (rightToggler.is(":visible")) {
        $("#map .ol-mouse-position")
                .css("margin-right", 50);
    } else {
        $("#map .ol-mouse-position")
                .css("margin-right", $(".sidebar-right").width());
    }
}
function isConstrained() {
    return $("div.mid").width() === $(window).width();
}
function applyInitialUIState() {
    if (isConstrained()) {
        $(".sidebar-left .sidebar-body").fadeOut('slide');
        $(".sidebar-right .sidebar-body").fadeOut('slide');
        $('.mini-submenu-left').fadeIn();
        $('.mini-submenu-right').fadeIn();
    }
}
$(function () {
    $('.sidebar-left .slide-submenu').on('click', function () {
        var thisEl = $(this);
        thisEl.closest('.sidebar-body').fadeOut('slide', function () {
            $('.mini-submenu-left').fadeIn();
            applyMargins();
        });
    });
    $('.mini-submenu-left').on('click', function () {
        var thisEl = $(this);
        $('.sidebar-left .sidebar-body').toggle('slide');
        thisEl.hide();
        applyMargins();
    });
    $('.sidebar-right .slide-submenu').on('click', function () {
        var thisEl = $(this);
        thisEl.closest('.sidebar-body').fadeOut('slide', function () {
            $('.mini-submenu-right').fadeIn();
            applyMargins();
        });
    });
    $('.mini-submenu-right').on('click', function () {
        var thisEl = $(this);
        $('.sidebar-right .sidebar-body').toggle('slide');
        thisEl.hide();
        applyMargins();
    });
    $(window).on("resize", applyMargins);
    applyInitialUIState();
    applyMargins();
});
