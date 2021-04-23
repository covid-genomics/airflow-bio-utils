from pathlib import Path

from poetry_publish.publish import poetry_publish

import airflow_bio_utils


def publish():
    poetry_publish(
        package_root=Path(airflow_bio_utils.__file__).parent.parent,
        version=airflow_bio_utils.__version__,
    )