#!/bin/bash

sed -i -e "s/http\:\/\/cdn\.mathjax/https\:\/\/cdn\.mathjax/g" `find ./ -name "*.html"`
