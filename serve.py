#!/usr/bin/env python3
"""Static server for report.html (chdirs itself so it works from anywhere)."""
import http.server
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
http.server.test(HandlerClass=http.server.SimpleHTTPRequestHandler, port=8090)
