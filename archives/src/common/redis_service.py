import asyncio
from collections import deque
import os, time, logging, json
import redis
from typing import Union
import os, pathlib, sys, argparse, threading, time
import subprocess
from multiprocessing.pool import ThreadPool

project_dir = str(pathlib.Path(__file__).resolve().parents[2])
sys.path.append(project_dir)

from src.common.logger import getLogger

logger = getLogger(__file__)


redis_host =   os.environ.get("REDIS_HOST", "localhost")  #"redis"
redis_user = os.environ.get("REDIS_USER", "dtcc_rmq" )
redis_password = os.environ.get("REDIS_PASSWORD", "dtcc_redis")
redis_port = os.environ.get("REDIS_PORT", "6879") # 6379

class RedisPubSub():
    def __init__(self,host:str=redis_host, port=redis_port) -> None:
        self.host = host
        self.port = int(port)
        self.client = self.get_connection(host=self.host, port=self.port)

    @staticmethod
    def get_connection(host="localhost", port=6379):
        return redis.Redis(host=host, port=port, password=redis_password, decode_responses=True)

    def test_redis(self) -> None:
        test_key = 'test'
        test_msg = 'It is working!'
        self.client.set(test_key, test_msg )
        
        assert test_msg == self.client.get(test_key)

    def publish(self, channel:str, message:str,retry_limit=5) -> bool:
        retries = 0
        published = False
        while True:
            try:
                rcvd = self.client.publish(channel, message)
                if rcvd >0:
                    logger.info("Success!")
                    published = True
                    break
            except redis.ConnectionError as e:
                logger.error(e)
                logger.error("Will attempt to retry")
                time.sleep(0.1)
            except Exception as e:
                logger.exception("from redis publish: "+str(e))
            
            if retries>=retry_limit:
                logger.error(f"Failed to send {message} to {channel} no receivers")
                break
            else:
                time.sleep(0.1)
            retries += 1
        return published       

    def subscribe_one(self,channel:str,wait_forever=True) -> Union[str, None]:
        message = None
        try:
            self.pubsub = self.client.pubsub(ignore_subscribe_messages=True)
            self.pubsub.subscribe(channel)
            retry_limit = 5
            retries = 0 
            while True:
                try:
                    message = self.pubsub.get_message()
                    if message and (not message['data'] == 1):
                        message = message['data']
                        break
                except redis.ConnectionError as e:
                    logger.error(e)
                    logger.error("Will attempt to retry")
                    
                except Exception as e:
                    logger.exception("from redis subscribe: "+str(e))
                
                if not wait_forever:
                    retries += 1
                    if retries>=retry_limit:
                        logger.error(f"Failed to receive message from {channel}")
                        break
        
                time.sleep(0.1)
        except redis.ConnectionError as e:
            logger.error(e)
        except Exception as e:
            logger.exception("from redis publish: "+ str(e))

        return message

    def subscribe(self,channel:str,callback=None) -> None:
        try:
            pubsub = self.client.pubsub(ignore_subscribe_messages=True)
            pubsub.subscribe(channel) 
            for raw_message in pubsub.listen():
                try:
                    if raw_message["type"] != "message":
                        continue
                    if raw_message and (not raw_message['data'] == 1):
                        message = raw_message["data"]
                        callback(message)
                except redis.ConnectionError as e:
                    logger.error(e)
                    logger.error("Will attempt to retry")
                    
                except Exception as e:
                    logger.exception("from redis subscribe: "+str(e))
                
        
                time.sleep(0.01)
        except redis.ConnectionError as e:
            logger.error(e)
        except Exception as e:
            logger.exception("from redis publish: "+ str(e))

    async def subscribe_async(self,channel:str):
        try:
            pubsub = self.client.pubsub(ignore_subscribe_messages=True)
            pubsub.subscribe(channel) 
            for raw_message in pubsub.listen():
                try:
                    if raw_message["type"] != "message":
                        continue
                    if raw_message and (not raw_message['data'] == 1):
                        message = raw_message["data"]
                        yield message
                except redis.ConnectionError as e:
                    logger.error(e)
                    logger.error("Will attempt to retry")
                    
                except Exception as e:
                    logger.exception("from redis subscribe: "+str(e))
                
        
                time.sleep(0.01)
        except redis.ConnectionError as e:
            logger.error(e)
        except Exception as e:
            logger.exception("from redis publish: "+ str(e))


def print_callback(msg:str):
    print(msg)


def test_send_receive(host, port):
    time_dt = subprocess.check_output(['date'])
    message=time_dt.decode('utf-8').replace('\n','')
    channel="core"

    rps = RedisPubSub(host=host,port=port)
    
    pool = ThreadPool(processes=1)

    async_result = pool.apply_async(rps.subscribe_one, args=(channel,))

    time.sleep(0.1)

    published = rps.publish(channel=channel,message=message)
    if published:
        print("published")

    incoming_msg = async_result.get()

    assert message == incoming_msg

    print(incoming_msg)

async def main():
    import asyncio

    import async_timeout

    import redis.asyncio as redis_async

    STOPWORD = "STOP"


    async def reader(channel: redis_async.client.PubSub):
        while True:
            try:
                async with async_timeout.timeout(1):
                    message = await channel.get_message(ignore_subscribe_messages=True)
                    if message is not None:
                        print(f"(Reader) Message Received: {message}")
                        if message["data"].decode() == STOPWORD:
                            print("(Reader) STOP")
                            break
                    await asyncio.sleep(0.1)
            except asyncio.TimeoutError:
                pass

    r = await redis_async.from_url(f"redis://:{redis_password}@localhost:6879/0")
    pubsub = r.pubsub()
    await pubsub.psubscribe("channel:1")

    

    await r.publish("channel:1", "Hello")
    await r.publish("channel:1", "World")
    await r.publish("channel:1", STOPWORD)

    future = asyncio.create_task(reader(pubsub))
    await future

    


if __name__=='__main__':
    asyncio.run(main())


    

if False:

    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers(dest='command')

    test = subparser.add_parser('test')
    test.add_argument("--host", "-H", type=str,default='localhost', help="hostname")
    test.add_argument("--port", "-P", type=int,default=6879, help="port")

    sub = subparser.add_parser('subscribe')
    sub.add_argument("--host", "-H", type=str,default='redis', help="hostname")
    sub.add_argument("--port", "-P", type=int,default=6379, help="port")
    sub.add_argument("--channel", "-C", type=str,default='dtcc-core', help="redis channel")


    pub = subparser.add_parser('publish')
    pub.add_argument("--host", "-H", type=str,default='redis', help="hostname")
    pub.add_argument("--port", "-P", type=int,default=6379, help="port")
    pub.add_argument("--channel", "-C", type=str,default='dtcc-core', help="redis channel")
    pub.add_argument("--message", "-m", type=str,default='', help="message to send")
    
    args = parser.parse_args()



    if args.command == 'test':
        test_send_receive(args.host, args.port)

    if args.command == 'subscribe':
        rps = RedisPubSub(host=args.host,port=args.port)
        rps.subscribe(channel=args.channel,callback=print_callback)

    if args.command == 'publish':
        rps = RedisPubSub(host=args.host,port=args.port)
        rps.publish(channel=args.channel,message=args.message, retry_limit=5)


    