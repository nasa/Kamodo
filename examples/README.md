# About

This directory contains example configuration files for kamodofied readers.

To test a given reader, use the kamodo command line tool:

```sh
kamodo config_override=/path/to/reader.yaml
```

For example, to test the HAPI reader, you would write:

```sh
kamodo config_override=hapi.yaml
```

!!! note
    Your python environment should all the dependencies for the given reader.

