import logging

from .readers import GSLC, RSLC

logging.basicConfig(
    format='{asctime} | {name:<15s} | {levelname:<7s} | {message}',
    datefmt='%Y-%m-%d %H:%M:%S',
    style='{',
    level=logging.INFO,
)