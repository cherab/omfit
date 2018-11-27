
from raysect.optical.observer import RGBPipeline2D, PowerPipeline2D, RadiancePipeline2D


def load_camera(config, world):

    if config['machine']['name'] == "MAST-U":

        from cherab.mastu.cameras import load_camera as load_mastu_camera

        camera = load_mastu_camera(config["observer"]["camera_id"], world,stride=config["observer"]["stride"])

    else:
        raise ValueError("Automatic camera loading is not supported for your machine.")

    pipelines = []

    if config["observer"]["display_progress"]:
        display_progress = True
    else:
        display_progress = False

    if config["observer"]["rgb_pipeline"]:
        rgb = RGBPipeline2D(display_unsaturated_fraction=0.96, name="sRGB",
                            display_progress=display_progress, display_update_time=15)
        pipelines.append(rgb)

    if config["observer"]["power_pipeline"]:
        power = PowerPipeline2D(display_unsaturated_fraction=0.96, name="Unfiltered Power (W)",
                                display_progress=display_progress, display_update_time=15)
        pipelines.append(power)

    if config["observer"]["radiance_pipeline"]:
        radiance = RadiancePipeline2D(display_unsaturated_fraction=0.96, name="Unfiltered Radiance (W/m^2/str)",
                                      display_progress=display_progress, display_update_time=15)
        pipelines.append(radiance)

    camera.pipelines = pipelines
    camera.pixel_samples = config["raytracing"]["pixel_samples"]

    return camera
